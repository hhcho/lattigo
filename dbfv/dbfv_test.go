package dbfv

import (
	"fmt"
	"log"
	"math/big"
	"testing"

	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"github.com/stretchr/testify/require"
)

func testString(opname string, parties uint64, params *bfv.Parameters) string {
	return fmt.Sprintf("%sparties=%d/LogN=%d/logQ=%d", opname, parties, params.LogN(), params.LogQP())
}

type dbfvTestContext struct {
	*dbfvContext

	params *bfv.Parameters

	prng utils.PRNG

	encoder bfv.Encoder
	kgen    *bfv.KeyGenerator

	sk0Shards []*bfv.SecretKey
	sk0       *bfv.SecretKey

	sk1       *bfv.SecretKey
	sk1Shards []*bfv.SecretKey

	pk0 *bfv.PublicKey
	pk1 *bfv.PublicKey

	encryptorPk0 bfv.Encryptor
	decryptorSk0 bfv.Decryptor
	decryptorSk1 bfv.Decryptor
	evaluator    bfv.Evaluator
}

type dbfvTestParameters struct {
	parties uint64

	contexts []*bfv.Parameters
}

var err error
var testParams = new(dbfvTestParameters)

func init() {
	testParams.parties = 3

	testParams.contexts = bfv.DefaultParams
}

func Test_DBFV(t *testing.T) {
	t.Run("PublicKeyGen", testPublicKeyGen)
	t.Run("RelinKeyGen", testRelinKeyGen)
	t.Run("RelinKeyGenNaive", testRelinKeyGenNaive)
	t.Run("KeySwitching", testKeyswitching)
	t.Run("PublicKeySwitching", testPublicKeySwitching)
	t.Run("RotKeyGenRotRows", testRotKeyGenRotRows)
	t.Run("RotKeyGenRotCols", testRotKeyGenRotCols)
	t.Run("Refresh", testRefresh)

}

func genDBFVTestContext(params *bfv.Parameters) (testCtx *dbfvTestContext) {

	testCtx = new(dbfvTestContext)

	testCtx.params = params
	testCtx.dbfvContext = newDbfvContext(params)

	testCtx.prng, err = utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	testCtx.encoder = bfv.NewEncoder(params)
	testCtx.evaluator = bfv.NewEvaluator(params)

	kgen := bfv.NewKeyGenerator(params)

	// SecretKeys
	testCtx.sk0Shards = make([]*bfv.SecretKey, testParams.parties)
	testCtx.sk1Shards = make([]*bfv.SecretKey, testParams.parties)
	tmp0 := testCtx.ringQP.NewPoly()
	tmp1 := testCtx.ringQP.NewPoly()

	for j := uint64(0); j < testParams.parties; j++ {
		testCtx.sk0Shards[j] = kgen.GenSecretKey()
		testCtx.sk1Shards[j] = kgen.GenSecretKey()
		testCtx.ringQP.Add(tmp0, testCtx.sk0Shards[j].Get(), tmp0)
		testCtx.ringQP.Add(tmp1, testCtx.sk1Shards[j].Get(), tmp1)
	}

	testCtx.sk0 = new(bfv.SecretKey)
	testCtx.sk1 = new(bfv.SecretKey)

	testCtx.sk0.Set(tmp0)
	testCtx.sk1.Set(tmp1)

	// Publickeys
	testCtx.pk0 = kgen.GenPublicKey(testCtx.sk0)
	testCtx.pk1 = kgen.GenPublicKey(testCtx.sk1)

	testCtx.encryptorPk0 = bfv.NewEncryptorFromPk(params, testCtx.pk0)
	testCtx.decryptorSk0 = bfv.NewDecryptor(params, testCtx.sk0)
	testCtx.decryptorSk1 = bfv.NewDecryptor(params, testCtx.sk1)

	return
}

func testPublicKeyGen(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards
		decryptorSk0 := testCtx.decryptorSk0
		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			crpGenerator := ring.NewUniformSampler(prng, testCtx.ringQP)
			crp := crpGenerator.ReadNew()

			type Party struct {
				*CKGProtocol
				s  *ring.Poly
				s1 CKGShare
			}

			ckgParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.CKGProtocol = NewCKGProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.s1 = p.AllocateShares()
				ckgParties[i] = p
			}
			P0 := ckgParties[0]

			// Each party creates a new CKGProtocol instance
			for i, p := range ckgParties {
				p.GenShare(p.s, crp, p.s1)
				if i > 0 {
					P0.AggregateShares(p.s1, P0.s1, P0.s1)
				}
			}

			pk := &bfv.PublicKey{}
			P0.GenPublicKey(P0.s1, crp, pk)

			// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
			encryptorTest := bfv.NewEncryptorFromPk(parameters, pk)

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorTest, t)

			verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)
		})
	}
}

func testRelinKeyGen(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk0 := testCtx.decryptorSk0
		evaluator := testCtx.evaluator

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*RKGProtocol
				u      *ring.Poly
				s      *ring.Poly
				share1 RKGShare
				share2 RKGShare
			}

			rkgParties := make([]*Party, parties)

			for i := range rkgParties {
				p := new(Party)
				p.RKGProtocol = NewEkgProtocol(parameters)
				p.u = p.RKGProtocol.NewEphemeralKey()
				p.s = sk0Shards[i].Get()
				p.share1, p.share2 = p.RKGProtocol.AllocateShares()
				rkgParties[i] = p
			}

			P0 := rkgParties[0]
			prng, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}

			crpGenerator := ring.NewUniformSampler(prng, testCtx.ringQP)
			crp := make([]*ring.Poly, parameters.Beta())

			for i := uint64(0); i < parameters.Beta(); i++ {
				crp[i] = crpGenerator.ReadNew()
			}

			// ROUND 1
			for i, p := range rkgParties {
				p.GenShareRoundOne(p.u, p.s, crp, p.share1)
				if i > 0 {
					P0.AggregateShareRoundOne(p.share1, P0.share1, P0.share1)
				}
			}

			//ROUND 2
			for i, p := range rkgParties {
				p.GenShareRoundTwo(P0.share1, p.u, p.s, crp, p.share2)
				if i > 0 {
					P0.AggregateShareRoundTwo(p.share2, P0.share2, P0.share2)
				}
			}

			evk := bfv.NewRelinKey(parameters, 1)
			P0.GenRelinearizationKey(P0.share1, P0.share2, evk)

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			for i := range coeffs {
				coeffs[i] *= coeffs[i]
				coeffs[i] %= testCtx.ringT.Modulus[0]
			}

			ciphertextMul := bfv.NewCiphertext(parameters, ciphertext.Degree()*2)
			evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

			res := bfv.NewCiphertext(parameters, 1)
			evaluator.Relinearize(ciphertextMul, evk, res)

			verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertextMul, t)
		})
	}
}

func testRelinKeyGenNaive(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		evaluator := testCtx.evaluator
		pk0 := testCtx.pk0
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk0 := testCtx.decryptorSk0
		sk0Shards := testCtx.sk0Shards

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*RKGProtocolNaive
				u      *ring.Poly
				s      *ring.Poly
				share1 RKGNaiveShareRoundOne
				share2 RKGNaiveShareRoundTwo
			}

			rkgParties := make([]*Party, parties)

			for i := range rkgParties {
				p := new(Party)
				p.RKGProtocolNaive = NewRKGProtocolNaive(parameters)
				p.s = sk0Shards[i].Get()
				p.share1, p.share2 = p.AllocateShares()
				rkgParties[i] = p
			}

			P0 := rkgParties[0]

			// ROUND 1
			for i, p := range rkgParties {
				rkgParties[i].GenShareRoundOne(p.s, pk0.Get(), p.share1)
				if i > 0 {
					P0.AggregateShareRoundOne(p.share1, P0.share1, P0.share1)
				}
			}

			// ROUND 2
			for i, p := range rkgParties {
				rkgParties[i].GenShareRoundTwo(P0.share1, p.s, pk0.Get(), p.share2)
				if i > 0 {
					P0.AggregateShareRoundTwo(p.share2, P0.share2, P0.share2)
				}
			}

			evk := bfv.NewRelinKey(parameters, 1)
			P0.GenRelinearizationKey(P0.share2, evk)

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			for i := range coeffs {
				coeffs[i] *= coeffs[i]
				coeffs[i] %= testCtx.ringT.Modulus[0]
			}

			ciphertextMul := bfv.NewCiphertext(parameters, ciphertext.Degree()*2)
			evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

			res := bfv.NewCiphertext(parameters, 1)
			evaluator.Relinearize(ciphertextMul, evk, res)

			verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertextMul, t)
		})
	}
}

func testKeyswitching(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards
		sk1Shards := testCtx.sk1Shards
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk1 := testCtx.decryptorSk1

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*CKSProtocol
				s0    *ring.Poly
				s1    *ring.Poly
				share CKSShare
			}

			cksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.CKSProtocol = NewCKSProtocol(parameters, 6.36)
				p.s0 = sk0Shards[i].Get()
				p.s1 = sk1Shards[i].Get()
				p.share = p.AllocateShare()
				cksParties[i] = p
			}
			P0 := cksParties[0]

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			// Each party creates its CKSProtocol instance with tmp = si-si'
			for i, p := range cksParties {
				p.GenShare(p.s0, p.s1, ciphertext, p.share)
				if i > 0 {
					P0.AggregateShares(p.share, P0.share, P0.share)
				}
			}

			ksCiphertext := bfv.NewCiphertext(parameters, 1)
			P0.KeySwitch(P0.share, ciphertext, ksCiphertext)

			verifyTestVectors(testCtx, decryptorSk1, coeffs, ksCiphertext, t)

			P0.KeySwitch(P0.share, ciphertext, ciphertext)

			verifyTestVectors(testCtx, decryptorSk1, coeffs, ciphertext, t)

		})
	}
}

func testPublicKeySwitching(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		sk0Shards := testCtx.sk0Shards
		pk1 := testCtx.pk1
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk1 := testCtx.decryptorSk1

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*PCKSProtocol
				s     *ring.Poly
				share PCKSShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.PCKSProtocol = NewPCKSProtocol(parameters, 6.36)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShares()
				pcksParties[i] = p
			}
			P0 := pcksParties[0]

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			ciphertextSwitched := bfv.NewCiphertext(parameters, 1)

			for i, p := range pcksParties {
				p.GenShare(p.s, pk1, ciphertext, p.share)
				if i > 0 {
					P0.AggregateShares(p.share, P0.share, P0.share)
				}
			}

			P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)

			verifyTestVectors(testCtx, decryptorSk1, coeffs, ciphertextSwitched, t)
		})
	}
}

func testRotKeyGenRotRows(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {

		testCtx := genDBFVTestContext(parameters)

		evaluator := testCtx.evaluator
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk0 := testCtx.decryptorSk0
		sk0Shards := testCtx.sk0Shards

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*RTGProtocol
				s     *ring.Poly
				share RTGShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RTGProtocol = NewRotKGProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShare()
				pcksParties[i] = p
			}
			P0 := pcksParties[0]
			prng, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}

			crpGenerator := ring.NewUniformSampler(prng, testCtx.ringQP)
			crp := make([]*ring.Poly, parameters.Beta())

			for i := uint64(0); i < parameters.Beta(); i++ {
				crp[i] = crpGenerator.ReadNew()
			}

			for i, p := range pcksParties {
				p.GenShare(bfv.RotationRow, 0, p.s, crp, &p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}

			rotkey := bfv.NewRotationKeys()
			P0.Finalize(P0.share, crp, rotkey)

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			evaluator.RotateRows(ciphertext, rotkey, ciphertext)

			coeffs = append(coeffs[testCtx.n>>1:], coeffs[:testCtx.n>>1]...)

			verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)

		})
	}

}

func testRotKeyGenRotCols(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		evaluator := testCtx.evaluator
		encryptorPk0 := testCtx.encryptorPk0
		decryptorSk0 := testCtx.decryptorSk0
		sk0Shards := testCtx.sk0Shards

		t.Run(testString("", parties, parameters), func(t *testing.T) {

			type Party struct {
				*RTGProtocol
				s     *ring.Poly
				share RTGShare
			}

			pcksParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RTGProtocol = NewRotKGProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShare()
				pcksParties[i] = p
			}

			P0 := pcksParties[0]
			prng, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}

			crpGenerator := ring.NewUniformSampler(prng, testCtx.ringQP)
			crp := make([]*ring.Poly, parameters.Beta())

			for i := uint64(0); i < parameters.Beta(); i++ {
				crp[i] = crpGenerator.ReadNew()
			}

			mask := (testCtx.n >> 1) - 1

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			receiver := bfv.NewCiphertext(parameters, ciphertext.Degree())

			for k := uint64(1); k < testCtx.n>>1; k <<= 1 {

				for i, p := range pcksParties {
					p.GenShare(bfv.RotationLeft, k, p.s, crp, &p.share)
					if i > 0 {
						P0.Aggregate(p.share, P0.share, P0.share)
					}
				}

				rotkey := bfv.NewRotationKeys()
				P0.Finalize(P0.share, crp, rotkey)

				evaluator.RotateColumns(ciphertext, k, rotkey, receiver)

				coeffsWant := make([]uint64, testCtx.n)

				for i := uint64(0); i < testCtx.n>>1; i++ {
					coeffsWant[i] = coeffs[(i+k)&mask]
					coeffsWant[i+(testCtx.n>>1)] = coeffs[((i+k)&mask)+(testCtx.n>>1)]
				}

				verifyTestVectors(testCtx, decryptorSk0, coeffsWant, receiver, t)
			}
		})
	}
}

func testRefresh(t *testing.T) {

	parties := testParams.parties

	for _, parameters := range testParams.contexts {
		testCtx := genDBFVTestContext(parameters)

		encryptorPk0 := testCtx.encryptorPk0
		sk0Shards := testCtx.sk0Shards
		encoder := testCtx.encoder
		decryptorSk0 := testCtx.decryptorSk0

		kgen := bfv.NewKeyGenerator(parameters)

		rlk := kgen.GenRelinKey(testCtx.sk0, 2)

		t.Run(fmt.Sprintf("N=%d/logQ=%d/Refresh", testCtx.n, testCtx.ringQP.ModulusBigint.BitLen()), func(t *testing.T) {

			type Party struct {
				*RefreshProtocol
				s       *ring.Poly
				share   RefreshShare
				ptShare *bfv.Plaintext
			}

			RefreshParties := make([]*Party, parties)
			for i := uint64(0); i < parties; i++ {
				p := new(Party)
				p.RefreshProtocol = NewRefreshProtocol(parameters)
				p.s = sk0Shards[i].Get()
				p.share = p.AllocateShares()
				p.ptShare = bfv.NewPlaintext(parameters)
				RefreshParties[i] = p
			}

			P0 := RefreshParties[0]
			prng, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}

			crpGenerator := ring.NewUniformSampler(prng, testCtx.ringQP)
			crp := crpGenerator.ReadNew()

			coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

			maxDepth := 0

			ciphertextTmp := ciphertext.CopyNew().Ciphertext()
			coeffsTmp := make([]uint64, len(coeffs))

			for i := range coeffs {
				coeffsTmp[i] = coeffs[i]
			}

			// Finds the maximum multiplicative depth
			for true {

				testCtx.evaluator.Relinearize(testCtx.evaluator.MulNew(ciphertextTmp, ciphertextTmp), rlk, ciphertextTmp)

				for j := range coeffsTmp {
					coeffsTmp[j] = ring.BRed(coeffsTmp[j], coeffsTmp[j], testCtx.ringT.Modulus[0], testCtx.ringT.BredParams[0])
				}

				if utils.EqualSliceUint64(coeffsTmp, encoder.DecodeUint(decryptorSk0.DecryptNew(ciphertextTmp))) {
					maxDepth++
				} else {
					break
				}
			}

			// Simulated added error of size Q/(T^2) and add it to the fresh ciphertext
			coeffsBigint := make([]*big.Int, testCtx.n)
			testCtx.ringQ.PolyToBigint(ciphertext.Value()[0], coeffsBigint)

			errorRange := new(big.Int).Set(testCtx.ringQ.ModulusBigint)
			errorRange.Quo(errorRange, testCtx.ringT.ModulusBigint)
			errorRange.Quo(errorRange, testCtx.ringT.ModulusBigint)

			for i := uint64(0); i < testCtx.n; i++ {
				coeffsBigint[i].Add(coeffsBigint[i], ring.RandInt(errorRange))
			}

			testCtx.ringQ.SetCoefficientsBigint(coeffsBigint, ciphertext.Value()[0])

			for i, p := range RefreshParties {
				p.GenShares(p.s, ciphertext, crp, p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}

			// We refresh the ciphertext with the simulated error
			P0.Decrypt(ciphertext, P0.share.RefreshShareDecrypt, P0.ptShare.Value()[0])      // Masked decryption
			P0.Recode(P0.ptShare.Value()[0], P0.ptShare.Value()[0])                          // Masked re-encoding
			P0.Recrypt(P0.ptShare.Value()[0], crp, P0.share.RefreshShareRecrypt, ciphertext) // Masked re-encryption$

			// The refresh also be called all at once with P0.Finalize(ciphertext, crp, P0.share, ctOut)

			// Square the refreshed ciphertext up to the maximum depth-1
			for i := 0; i < maxDepth-1; i++ {

				testCtx.evaluator.Relinearize(testCtx.evaluator.MulNew(ciphertext, ciphertext), rlk, ciphertext)

				for j := range coeffs {
					coeffs[j] = ring.BRed(coeffs[j], coeffs[j], testCtx.ringT.Modulus[0], testCtx.ringT.BredParams[0])
				}
			}

			//Decrypts and compare
			require.True(t, utils.EqualSliceUint64(coeffs, encoder.DecodeUint(decryptorSk0.DecryptNew(ciphertext))))
		})
	}
}

func newTestVectors(contextParams *dbfvTestContext, encryptor bfv.Encryptor, t *testing.T) (coeffs []uint64, plaintext *bfv.Plaintext, ciphertext *bfv.Ciphertext) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	uniformSampler := ring.NewUniformSampler(prng, contextParams.ringT)

	coeffsPol := uniformSampler.ReadNew()
	plaintext = bfv.NewPlaintext(contextParams.params)
	contextParams.encoder.EncodeUint(coeffsPol.Coeffs[0], plaintext)
	ciphertext = encryptor.EncryptNew(plaintext)
	return coeffsPol.Coeffs[0], plaintext, ciphertext
}

func verifyTestVectors(contextParams *dbfvTestContext, decryptor bfv.Decryptor, coeffs []uint64, ciphertext *bfv.Ciphertext, t *testing.T) {
	require.True(t, utils.EqualSliceUint64(coeffs, contextParams.encoder.DecodeUint(decryptor.DecryptNew(ciphertext))))
}

func Test_Marshalling(t *testing.T) {
	params := bfv.DefaultParams[bfv.PN14QP438]

	//verify if the un.marshalling works properly
	dbfvCtx := newDbfvContext(params)
	KeyGenerator := bfv.NewKeyGenerator(params)
	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	if err != nil {
		panic(err)
	}
	crsGen := ring.NewUniformSampler(prng, dbfvCtx.ringQP)
	sk := KeyGenerator.GenSecretKey()
	crs := crsGen.ReadNew()
	ringQ := dbfvCtx.ringQ
	contextPKeys := dbfvCtx.ringP

	Ciphertext := bfv.NewCiphertextRandom(prng, params, 1)

	t.Run(fmt.Sprintf("CPK/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {
		keygenProtocol := NewCKGProtocol(params)
		KeyGenShareBefore := keygenProtocol.AllocateShares()
		keygenProtocol.GenShare(sk.Get(), crs, KeyGenShareBefore)
		//now we marshall it
		data, err := KeyGenShareBefore.MarshalBinary()

		if err != nil {
			log.Fatal("Could not marshal the CKGShare : ", err)

		}

		KeyGenShareAfter := new(CKGShare)
		err = KeyGenShareAfter.UnmarshalBinary(data)
		if err != nil {
			log.Fatal("Could not unmarshal the CKGShare : ", err)

		}

		//comparing the results
		require.Equal(t, KeyGenShareBefore.GetDegree(), KeyGenShareAfter.GetDegree())
		require.Equal(t, KeyGenShareBefore.GetLenModuli(), KeyGenShareAfter.GetLenModuli())

		moduli := KeyGenShareBefore.GetLenModuli()
		require.Equal(t, KeyGenShareAfter.Coeffs[:moduli], KeyGenShareBefore.Coeffs[:moduli])
	})

	t.Run(fmt.Sprintf("PCKS/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {
		//Check marshalling for the PCKS

		KeySwitchProtocol := NewPCKSProtocol(params, dbfvCtx.params.Sigma())
		SwitchShare := KeySwitchProtocol.AllocateShares()
		pk := KeyGenerator.GenPublicKey(sk)
		KeySwitchProtocol.GenShare(sk.Get(), pk, Ciphertext, SwitchShare)

		data, err := SwitchShare.MarshalBinary()
		require.NoError(t, err)

		SwitchShareReceiver := new(PCKSShare)
		err = SwitchShareReceiver.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := 0; i < 2; i++ {
			//compare the shares.
			ringBefore := SwitchShare[i]
			ringAfter := SwitchShareReceiver[i]
			require.Equal(t, ringBefore.GetDegree(), ringAfter.GetDegree())
			moduli := ringAfter.GetLenModuli()
			require.Equal(t, ringAfter.Coeffs[:moduli], ringBefore.Coeffs[:moduli])
		}
	})

	t.Run(fmt.Sprintf("CKS/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {

		//Now for CKSShare ~ its similar to PKSShare
		cksp := NewCKSProtocol(params, dbfvCtx.params.Sigma())
		cksshare := cksp.AllocateShare()
		skIn := KeyGenerator.GenSecretKey()
		skOut := KeyGenerator.GenSecretKey()
		cksp.GenShare(skIn.Get(), skOut.Get(), Ciphertext, cksshare)

		data, err := cksshare.MarshalBinary()
		require.NoError(t, err)
		cksshareAfter := new(CKSShare)
		err = cksshareAfter.UnmarshalBinary(data)
		require.NoError(t, err)

		//now compare both shares.

		require.Equal(t, cksshare.GetDegree(), cksshareAfter.GetDegree())
		require.Equal(t, cksshare.GetLenModuli(), cksshareAfter.GetLenModuli())

		moduli := cksshare.GetLenModuli()
		require.Equal(t, cksshare.Coeffs[:moduli], cksshareAfter.Coeffs[:moduli])
	})

	t.Run(fmt.Sprintf("Refresh/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {

		//testing refresh shares
		refreshproto := NewRefreshProtocol(params)
		refreshshare := refreshproto.AllocateShares()
		refreshproto.GenShares(sk.Get(), Ciphertext, crs, refreshshare)

		data, err := refreshshare.MarshalBinary()
		if err != nil {
			log.Fatal("Could not marshal RefreshShare", err)
		}
		resRefreshShare := new(RefreshShare)
		err = resRefreshShare.UnmarshalBinary(data)

		if err != nil {
			log.Fatal("Could not unmarshal RefreshShare", err)
		}
		for i, r := range refreshshare.RefreshShareDecrypt.Coeffs {
			if !utils.EqualSliceUint64(resRefreshShare.RefreshShareDecrypt.Coeffs[i], r) {
				log.Fatal("Resulting of marshalling not the same as original : RefreshShare")
			}

		}
		for i, r := range refreshshare.RefreshShareRecrypt.Coeffs {
			if !utils.EqualSliceUint64(resRefreshShare.RefreshShareRecrypt.Coeffs[i], r) {
				log.Fatal("Resulting of marshalling not the same as original : RefreshShare")
			}

		}
	})

	t.Run(fmt.Sprintf("RTG/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {

		//check RTGShare
		prng, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}
		crpGenerator := ring.NewUniformSampler(prng, dbfvCtx.ringQP)
		modulus := (dbfvCtx.ringQ.Modulus)
		crp := make([]*ring.Poly, len(modulus))
		for j := 0; j < len(modulus); j++ {
			crp[j] = crpGenerator.ReadNew() //make([]*ring.Poly, bitLog)

		}

		rotProto := NewRotKGProtocol(params)
		rtgShare := rotProto.AllocateShare()
		rotProto.GenShare(1, 64, sk.Get(), crp, &rtgShare)

		data, err := rtgShare.MarshalBinary()
		if err != nil {
			log.Fatal("could not marshal RTGshare :", err)
		}

		resRTGShare := new(RTGShare)
		err = resRTGShare.UnmarshalBinary(data)
		if err != nil {
			log.Fatal("Could not unmarshal RTGShare: ", err)
		}

		if resRTGShare.Type != rtgShare.Type || resRTGShare.K != rtgShare.K || len(resRTGShare.Value) != len(rtgShare.Value) {
			log.Fatal("result after marshalling is not the same as before marshalling for RTGSahre")
		}

		for i, val := range rtgShare.Value {
			ring1 := val
			ring2 := resRTGShare.Value[i]
			if len(ring1.Coeffs) != len(ring2.Coeffs) {
				log.Fatal("result after marshalling is not the same as before marshalling for RTGSahre")
			}

			for j, elem := range ring1.Coeffs {
				if !utils.EqualSliceUint64(ring2.Coeffs[j], elem) {
					log.Fatal("result after marshalling is not the same as before marshalling for RTGSahre")

				}
			}

		}
	})

}

func Test_Relin_Marshalling(t *testing.T) {
	params := bfv.DefaultParams[bfv.PN14QP438]

	dbfvCtx := newDbfvContext(params)
	ringQ := dbfvCtx.ringQ
	contextPKeys := dbfvCtx.ringP
	modulus := dbfvCtx.ringQ.Modulus
	prng, err := utils.NewKeyedPRNG(nil)
	if err != nil {
		panic(err)
	}

	crpGenerator := ring.NewUniformSampler(prng, dbfvCtx.ringQP)

	crp := make([]*ring.Poly, len(modulus))
	for j := 0; j < len(modulus); j++ {
		crp[j] = crpGenerator.ReadNew() //make([]*ring.Poly, bitLog)
		//for u := uint64(0); u < bitLog; u++ {
		//	crp[j][u] = crpGenerator.ClockUniformNew()
		//}
	}

	t.Run(fmt.Sprintf("RLKG/N=%d/limbQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(contextPKeys.Modulus)), func(t *testing.T) {

		rlk := NewEkgProtocol(params)
		u := rlk.NewEphemeralKey()
		sk := bfv.NewKeyGenerator(params).GenSecretKey()

		r1, r2 := rlk.AllocateShares()
		rlk.GenShareRoundOne(u, sk.Get(), crp, r1)
		data, err := r1.MarshalBinary()
		require.NoError(t, err)

		r1After := new(RKGShare)
		err = r1After.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := 0; i < (len(r1)); i++ {
			a := r1[i][0]
			b := (*r1After)[i][0]
			moduli := a.GetLenModuli()
			require.Equal(t, a.Coeffs[:moduli], b.Coeffs[:moduli])
		}

		rlk.GenShareRoundTwo(r1, u, sk.Get(), crp, r2)

		data, err = r2.MarshalBinary()
		require.NoError(t, err)

		r2After := new(RKGShare)
		err = r2After.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := 0; i < (len(r2)); i++ {
			for idx := 0; idx < 2; idx++ {
				a := r2[i][idx]
				b := (*r2After)[i][idx]
				moduli := a.GetLenModuli()
				require.Equal(t, a.Coeffs[:moduli], b.Coeffs[:moduli])
			}

		}
	})

}

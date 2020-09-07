// Package ring implelents a RNS-accelerated modular arithmetic operations for polynomials, including: RNS basis extension; RNS rescaling;  number theoretic transform (NTT); uniform, Gaussian and ternary sampling.
package ring

import (
	"bytes"
	"encoding/gob"
	"errors"
	"math/big"
	"math/bits"

	"github.com/ldsec/lattigo/utils"
)

// Ring serves as an evaluation context for ring elements.
// It is a structure holding all parameters and pre-computations that are necessary for the arithmetic operations in the defined RNS polynomial quotient ring.
type Ring struct {

	// Polynomial nb.Coefficients
	N uint64

	// Moduli
	Modulus []uint64

	// 2^bit_length(Qi) - 1
	Mask []uint64

	// Determines if NTT can be used with the current context.
	AllowsNTT bool

	// Product of the Moduli
	ModulusBigint *big.Int

	// Fast reduction parameters
	BredParams [][]uint64
	MredParams []uint64

	RescaleParams [][]uint64

	//NTT Parameters
	PsiMont    []uint64 //2nth primitive root in Montgomery form
	PsiInvMont []uint64 //2nth inverse primitive root in Montgomery form

	NttPsi    [][]uint64 //powers of the inverse of the 2nth primitive root in Montgomery form (in bitreversed order)
	NttPsiInv [][]uint64 //powers of the inverse of the 2nth primitive root in Montgomery form (in bitreversed order)
	NttNInv   []uint64   //[N^-1] mod Qi in Montgomery form
}

// NewRing creates a new ringContex with the given parameters. Checks that N is a power of 2 and that the moduli are NTT compliant.
func NewRing(N uint64, Moduli []uint64) (r *Ring, err error) {
	r = new(Ring)
	r.setParameters(N, Moduli)
	return r, r.genNTTParams()
}

// setParameters initialises a *Ring by setting the required pre-computed values (except for the NTT-related values which are set by the
// genNTTParams function).
func (r *Ring) setParameters(N uint64, Modulus []uint64) {

	// Checks if N is a power of 2
	if (N&(N-1)) != 0 && N != 0 {
		panic("invalid ring degree (must be a power of 2)")
	}

	r.AllowsNTT = false

	r.N = N

	r.Modulus = make([]uint64, len(Modulus))
	r.Mask = make([]uint64, len(Modulus))

	for i, qi := range Modulus {
		r.Modulus[i] = qi
		r.Mask[i] = (1 << uint64(bits.Len64(qi))) - 1
	}

	//Computes the bigQ
	r.ModulusBigint = NewInt(1)
	for _, qi := range r.Modulus {
		r.ModulusBigint.Mul(r.ModulusBigint, NewUint(qi))
	}

	// Computes the fast reduction parameters
	r.BredParams = make([][]uint64, len(r.Modulus))
	r.MredParams = make([]uint64, len(r.Modulus))

	for i, qi := range r.Modulus {

		//Computes the fast modular reduction parameters for the Context
		r.BredParams[i] = BRedParams(qi)

		// If qi is not a power of 2, we can compute the MRedParams (else it should not
		// because it will return an error and there is no valid Montgomery form mod a power of 2)
		if (qi&(qi-1)) != 0 && qi != 0 {
			r.MredParams[i] = MRedParams(qi)
		}
	}

}

// genNTTParams checks that N has been correctly initialized, and checks that each moduli is a prime congruent to 1 mod 2N (i.e. allowing NTT).
// Then it computes the variables required for the NTT. ValidateParameters purpose is to validate that the moduli allow the NTT and compute the
// NTT parameters.
func (r *Ring) genNTTParams() error {

	if r.AllowsNTT {
		return nil
	}

	if r.N == 0 || r.Modulus == nil {
		panic("error : invalid context parameters (missing)")
	}

	// Checks if each qi is Prime and if qi = 1 mod 2n
	for _, qi := range r.Modulus {
		if IsPrime(qi) == false || qi&((r.N<<1)-1) != 1 {
			r.AllowsNTT = false
			return errors.New("warning : provided modulus does not allow NTT")
		}
	}

	r.RescaleParams = make([][]uint64, len(r.Modulus)-1)

	for j := len(r.Modulus) - 1; j > 0; j-- {

		r.RescaleParams[j-1] = make([]uint64, j)

		for i := 0; i < j; i++ {

			r.RescaleParams[j-1][i] = MForm(ModExp(r.Modulus[j], r.Modulus[i]-2, r.Modulus[i]), r.Modulus[i], r.BredParams[i])
		}
	}

	r.PsiMont = make([]uint64, len(r.Modulus))
	r.PsiInvMont = make([]uint64, len(r.Modulus))
	r.NttPsi = make([][]uint64, len(r.Modulus))
	r.NttPsiInv = make([][]uint64, len(r.Modulus))
	r.NttNInv = make([]uint64, len(r.Modulus))

	bitLenofN := uint64(bits.Len64(r.N) - 1)

	for i, qi := range r.Modulus {

		//2.1 Computes N^(-1) mod Q in Montgomery form
		r.NttNInv[i] = MForm(ModExp(r.N, qi-2, qi), qi, r.BredParams[i])

		//2.2 Computes Psi and PsiInv in Montgomery form
		r.NttPsi[i] = make([]uint64, r.N)
		r.NttPsiInv[i] = make([]uint64, r.N)

		//Finds a 2nth primitive Root
		g := primitiveRoot(qi)

		_2n := uint64(r.N << 1)

		power := (qi - 1) / _2n
		powerInv := (qi - 1) - power

		//Computes Psi and PsiInv in Montgomery Form
		PsiMont := MForm(ModExp(g, power, qi), qi, r.BredParams[i])
		PsiInvMont := MForm(ModExp(g, powerInv, qi), qi, r.BredParams[i])

		r.PsiMont[i] = PsiMont
		r.PsiInvMont[i] = PsiInvMont

		r.NttPsi[i][0] = MForm(1, qi, r.BredParams[i])
		r.NttPsiInv[i][0] = MForm(1, qi, r.BredParams[i])

		// Computes nttPsi[j] = nttPsi[j-1]*Psi and nttPsiInv[j] = nttPsiInv[j-1]*PsiInv
		for j := uint64(1); j < r.N; j++ {

			indexReversePrev := utils.BitReverse64(j-1, bitLenofN)
			indexReverseNext := utils.BitReverse64(j, bitLenofN)

			r.NttPsi[i][indexReverseNext] = MRed(r.NttPsi[i][indexReversePrev], PsiMont, qi, r.MredParams[i])
			r.NttPsiInv[i][indexReverseNext] = MRed(r.NttPsiInv[i][indexReversePrev], PsiInvMont, qi, r.MredParams[i])
		}
	}

	r.AllowsNTT = true

	return nil
}

// Used to export the context. Minimal information to recover the full context.
type smallContext struct {
	N       uint64
	Modulus []uint64
}

// MarshalBinary encodes the target ring context on a slice of bytes.
func (r *Ring) MarshalBinary() ([]byte, error) {

	parameters := smallContext{r.N, r.Modulus}

	var buf bytes.Buffer
	enc := gob.NewEncoder(&buf)
	if err := enc.Encode(parameters); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

// UnmarshalBinary decodes slice of bytes on the target ring context.
func (r *Ring) UnmarshalBinary(data []byte) error {

	parameters := smallContext{}

	reader := bytes.NewReader(data)
	dec := gob.NewDecoder(reader)
	if err := dec.Decode(&parameters); err != nil {
		return err
	}

	r.setParameters(parameters.N, parameters.Modulus)
	r.genNTTParams()

	return nil
}

// NewPoly create a new polynomial with all coefficients set to 0.
func (r *Ring) NewPoly() *Poly {
	p := new(Poly)

	p.Coeffs = make([][]uint64, len(r.Modulus))
	for i := 0; i < len(r.Modulus); i++ {
		p.Coeffs[i] = make([]uint64, r.N)
	}

	return p
}

// NewPolyLvl create a new polynomial with all coefficients set to 0.
func (r *Ring) NewPolyLvl(level uint64) *Poly {
	p := new(Poly)

	p.Coeffs = make([][]uint64, level+1)
	for i := uint64(0); i < level+1; i++ {
		p.Coeffs[i] = make([]uint64, r.N)
	}

	return p
}

// SetCoefficientsInt64 sets the coefficients of p1 from an int64 array.
func (r *Ring) SetCoefficientsInt64(coeffs []int64, p1 *Poly) {
	for i, coeff := range coeffs {
		for j, Qi := range r.Modulus {
			p1.Coeffs[j][i] = CRed(uint64((coeff%int64(Qi) + int64(Qi))), Qi)
		}
	}
}

// SetCoefficientsUint64 sets the coefficients of p1 from an uint64 array.
func (r *Ring) SetCoefficientsUint64(coeffs []uint64, p1 *Poly) {
	for i, coeff := range coeffs {
		for j, Qi := range r.Modulus {
			p1.Coeffs[j][i] = coeff % Qi
		}
	}
}

// SetCoefficientsString parses an array of string as Int variables, and sets the
// coefficients of p1 with this Int variables.
func (r *Ring) SetCoefficientsString(coeffs []string, p1 *Poly) {
	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i, Qi := range r.Modulus {
		QiBigint.SetUint64(Qi)
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(NewIntFromString(coeff), QiBigint).Uint64()
		}
	}
}

// SetCoefficientsBigint sets the coefficients of p1 from an array of Int variables.
func (r *Ring) SetCoefficientsBigint(coeffs []*big.Int, p1 *Poly) {
	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i, Qi := range r.Modulus {
		QiBigint.SetUint64(Qi)
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(coeff, QiBigint).Uint64()

		}
	}
}

// SetCoefficientsBigintLvl sets the coefficients of p1 from an array of Int variables.
func (r *Ring) SetCoefficientsBigintLvl(level uint64, coeffs []*big.Int, p1 *Poly) {

	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i := uint64(0); i < level+1; i++ {
		QiBigint.SetUint64(r.Modulus[i])
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(coeff, QiBigint).Uint64()

		}
	}
}

// PolyToString reconstructs p1 and returns the result in an array of string.
func (r *Ring) PolyToString(p1 *Poly) []string {

	coeffsBigint := make([]*big.Int, r.N)
	r.PolyToBigint(p1, coeffsBigint)
	coeffsString := make([]string, len(coeffsBigint))

	for i := range coeffsBigint {
		coeffsString[i] = coeffsBigint[i].String()
	}

	return coeffsString
}

// PolyToBigint reconstructs p1 and returns the result in an array of Int.
func (r *Ring) PolyToBigint(p1 *Poly, coeffsBigint []*big.Int) {

	var qi, level uint64

	level = uint64(len(p1.Coeffs) - 1)

	crtReconstruction := make([]*big.Int, level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := NewUint(1)

	for i := uint64(0); i < level+1; i++ {

		qi = r.Modulus[i]
		QiB.SetUint64(qi)

		modulusBigint.Mul(modulusBigint, QiB)

		crtReconstruction[i] = new(big.Int)
		crtReconstruction[i].Quo(r.ModulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	for x := uint64(0); x < r.N; x++ {

		tmp.SetUint64(0)
		coeffsBigint[x] = new(big.Int)

		for i := uint64(0); i < level+1; i++ {
			coeffsBigint[x].Add(coeffsBigint[x], tmp.Mul(NewUint(p1.Coeffs[i][x]), crtReconstruction[i]))
		}

		coeffsBigint[x].Mod(coeffsBigint[x], modulusBigint)
	}
}

// PolyToBigint reconstructs p1 and returns the result in an pre-allocated array of Int.
func (r *Ring) PolyToBigintNoAlloc(p1 *Poly, coeffsBigint []*big.Int) {

	var qi, level uint64

	level = uint64(len(p1.Coeffs) - 1)

	crtReconstruction := make([]*big.Int, level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := NewUint(1)

	for i := uint64(0); i < level+1; i++ {

		qi = r.Modulus[i]
		QiB.SetUint64(qi)

		modulusBigint.Mul(modulusBigint, QiB)

		crtReconstruction[i] = new(big.Int)
		crtReconstruction[i].Quo(r.ModulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	for x := uint64(0); x < r.N; x++ {

		tmp.SetUint64(0)

		for i := uint64(0); i < level+1; i++ {
			coeffsBigint[x].Add(coeffsBigint[x], tmp.Mul(NewUint(p1.Coeffs[i][x]), crtReconstruction[i]))
		}

		coeffsBigint[x].Mod(coeffsBigint[x], modulusBigint)
	}
}

// Equal checks if p1 = p2 in the given context.
func (r *Ring) Equal(p1, p2 *Poly) bool {

	for i := 0; i < len(r.Modulus); i++ {
		if len(p1.Coeffs[i]) != len(p2.Coeffs[i]) {
			return false
		}
	}

	r.Reduce(p1, p1)
	r.Reduce(p2, p2)

	for i := 0; i < len(r.Modulus); i++ {
		for j := uint64(0); j < r.N; j++ {
			if p1.Coeffs[i][j] != p2.Coeffs[i][j] {
				return false
			}
		}
	}

	return true
}

// EqualLvl checks if p1 = p2 in the given context.
func (r *Ring) EqualLvl(level uint64, p1, p2 *Poly) bool {

	for i := uint64(0); i < level+1; i++ {
		if len(p1.Coeffs[i]) != len(p2.Coeffs[i]) {
			return false
		}
	}

	r.ReduceLvl(level, p1, p1)
	r.ReduceLvl(level, p2, p2)

	for i := uint64(0); i < level+1; i++ {
		for j := uint64(0); j < r.N; j++ {
			if p1.Coeffs[i][j] != p2.Coeffs[i][j] {
				return false
			}
		}
	}

	return true
}

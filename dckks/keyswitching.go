package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {
	dckksContext *Context

	sigmaSmudging         float64
	gaussianSamplerSmudge *ring.KYSampler

	tmp      *ring.Poly
	tmpDelta *ring.Poly
	hP       *ring.Poly

	baseconverter *ring.FastBasisExtender
}

// CKSShare is a struct holding a share of the CKS protocol.
type CKSShare *ring.Poly

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params *ckks.Parameters, sigmaSmudging float64) (cks *CKSProtocol) {

	if !params.IsValid() {
		panic("cannot NewCKSProtocol : params not valid (check if they where generated properly)")
	}

	cks = new(CKSProtocol)

	dckksContext := NewContext(params)

	cks.dckksContext = dckksContext

	cks.gaussianSamplerSmudge = dckksContext.ContextQP.NewKYSampler(sigmaSmudging, int(6*sigmaSmudging))

	cks.tmp = dckksContext.ContextQP.NewPoly()
	cks.tmpDelta = dckksContext.ContextQ.NewPoly()
	cks.hP = dckksContext.ContextP.NewPoly()

	cks.baseconverter = ring.NewFastBasisExtender(dckksContext.ContextQ, dckksContext.ContextP)

	return cks
}

// AllocateShare allocates the share of the CKS protocol.
func (cks *CKSProtocol) AllocateShare() CKSShare {
	return cks.dckksContext.ContextQ.NewPoly()
}

// GenShare is the first and unique round of the CKSProtocol protocol. Each party holding a ciphertext ctx encrypted under a collective publick-key must
// compute the following :
//
// [(skInput_i - skOutput_i) * ctx[0] + e_i]
//
// Each party then broadcasts the result of this computation to the other j-1 parties.
func (cks *CKSProtocol) GenShare(skInput, skOutput *ring.Poly, ct *ckks.Ciphertext, shareOut CKSShare) {

	cks.dckksContext.ContextQ.Sub(skInput, skOutput, cks.tmpDelta)

	cks.genShareDelta(cks.tmpDelta, ct, shareOut)
}

func (cks *CKSProtocol) genShareDelta(skDelta *ring.Poly, ct *ckks.Ciphertext, shareOut CKSShare) {

	contextQ := cks.dckksContext.ContextQ
	contextP := cks.dckksContext.ContextP

	contextQ.MulCoeffsMontgomeryLvl(ct.Level(), ct.Value()[1], skDelta, shareOut)

	contextQ.MulScalarBigintLvl(ct.Level(), shareOut, contextP.ModulusBigint, shareOut)

	// TODO : improve by only computing the NTT for the required primes
	cks.gaussianSamplerSmudge.SampleNTT(cks.tmp)
	contextQ.AddLvl(ct.Level(), shareOut, cks.tmp, shareOut)

	for x, i := 0, uint64(len(contextQ.Modulus)); i < uint64(len(cks.dckksContext.ContextQP.Modulus)); x, i = x+1, i+1 {
		tmp0 := cks.tmp.Coeffs[i]
		tmp1 := cks.hP.Coeffs[x]
		for j := uint64(0); j < contextQ.N; j++ {
			tmp1[j] += tmp0[j]
		}
	}

	cks.baseconverter.ModDownSplitedNTTPQ(ct.Level(), shareOut, cks.hP, shareOut)

	cks.hP.Zero()
	cks.tmp.Zero()
}

// AggregateShares is the second part of the unique round of the CKSProtocol protocol. Upon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *CKSProtocol) AggregateShares(share1, share2, shareOut CKSShare) {
	cks.dckksContext.ContextQ.AddLvl(uint64(len(share1.Coeffs)-1), share1, share2, shareOut)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined CKSShare, ct *ckks.Ciphertext, ctOut *ckks.Ciphertext) {
	ctOut.SetScale(ct.Scale())
	cks.dckksContext.ContextQ.AddLvl(ct.Level(), ct.Value()[0], combined, ctOut.Value()[0])
	cks.dckksContext.ContextQ.CopyLvl(ct.Level(), ct.Value()[1], ctOut.Value()[1])
}

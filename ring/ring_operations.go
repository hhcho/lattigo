package ring

import (
	"math/big"
	"math/bits"

	"github.com/ldsec/lattigo/utils"
)

// Add adds p1 to p2 coefficient wise and writes the result on p3.
func (r *Ring) Add(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = CRed(p1tmp[j]+p2tmp[j], qi)
		}
	}
}

// AddLvl adds p1 to p2 coefficient wise for the moduli
// q_0 up to q_level and writes the result on p3.
func (r *Ring) AddLvl(level uint64, p1, p2, p3 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = CRed(p1tmp[j]+p2tmp[j], qi)
		}
	}
}

// AddNoMod adds p1 to p2 coefficient wise without
// modular reduction and writes the result on p3.
func (r *Ring) AddNoMod(p1, p2, p3 *Poly) {
	for i := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = p1tmp[j] + p2tmp[j]
		}
	}
}

// AddNoModLvl adds p1 to p2 coefficient wise without modular reduction
// for the moduli q_0 up to q_level and writes the result on p3.
func (r *Ring) AddNoModLvl(level uint64, p1, p2, p3 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = p1tmp[j] + p2tmp[j]
		}
	}
}

// Sub subtracts p2 to p1 coefficient wise and writes the result on p3.
func (r *Ring) Sub(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = CRed((p1tmp[j]+qi)-p2tmp[j], qi)
		}
	}
}

// SubLvl subtracts p2 to p1 coefficient wise and writes the result on p3.
func (r *Ring) SubLvl(level uint64, p1, p2, p3 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = CRed((p1tmp[j]+qi)-p2tmp[j], qi)
		}
	}
}

// SubNoMod subtracts p2 to p1 coefficient wise without
// modular reduction and returns the result on p3.
func (r *Ring) SubNoMod(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = (p1tmp[j] + qi) - p2tmp[j]
		}
	}
}

// SubNoModLvl subtracts p2 to p1 coefficient wise without modular reduction
// for the moduli q_0 up to q_level and writes the result on p3.
func (r *Ring) SubNoModLvl(level uint64, p1, p2, p3 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = (p1tmp[j] + qi) - p2tmp[j]
		}
	}
}

// Neg sets all coefficients of p1 to there additive inverse and writes the result on p2.
func (r *Ring) Neg(p1, p2 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = qi - p1tmp[j]
		}
	}
}

// NegLvl sets the coefficients of p1 to there additive inverse for
// the moduli q_0 up to q_level and writes the result on p2.
func (r *Ring) NegLvl(level uint64, p1, p2 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = qi - p1tmp[j]
		}
	}
}

// Reduce applies a modular reduction on the coefficients of p1 and writes the result on p2.
func (r *Ring) Reduce(p1, p2 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = BRedAdd(p1tmp[j], qi, bredParams)
		}
	}
}

// ReduceLvl applies a modular reduction over the coefficients of p1
// for the moduli q_0 up to q_level and writes the result on p2.
func (r *Ring) ReduceLvl(level uint64, p1, p2 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = BRedAdd(p1tmp[j], qi, bredParams)
		}
	}
}

// Mod applies a modular reduction by m over the coefficients of p1 and writes the result on p2.
func (r *Ring) Mod(p1 *Poly, m uint64, p2 *Poly) {
	params := BRedParams(m)
	for i := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = BRedAdd(p1tmp[j], m, params)
		}
	}
}

// AND applies a logical AND of m to the coefficients of p1 and writes the result on p2.
func (r *Ring) AND(p1 *Poly, m uint64, p2 *Poly) {
	for i := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = p1tmp[j] & m
		}
	}
}

// OR applies a logical OR of m to the coefficients of p1 and writes the result on p2.
func (r *Ring) OR(p1 *Poly, m uint64, p2 *Poly) {
	for i := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = p1tmp[j] | m
		}
	}
}

// XOR applies a logical XOR of m to the coefficients of p1 and writes the result on p2.
func (r *Ring) XOR(p1 *Poly, m uint64, p2 *Poly) {
	for i := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = p1tmp[j] ^ m
		}
	}
}

// MulCoeffs multiplies p1 by p2 coefficient wise and uses a
// Barrett modular reduction and writes the result on p3.
func (r *Ring) MulCoeffs(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = BRed(p1tmp[j], p2tmp[j], qi, bredParams)
		}
	}
}

// MulCoeffsAndAdd multiplies p1 by p2 coefficient wise with
// a Barret modular reduction and adds the result to p3.
func (r *Ring) MulCoeffsAndAdd(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = CRed(p3tmp[j]+BRed(p1tmp[j], p2tmp[j], qi, bredParams), qi)
		}
	}
}

// MulCoeffsAndAddNoMod multiplies p1 by p2 coefficient wise with a Barrett
// modular reduction and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsAndAddNoMod(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] += BRed(p1tmp[j], p2tmp[j], qi, bredParams)
		}
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient wise with a
// Montgomery modular reduction and returns the result on p3.
func (r *Ring) MulCoeffsMontgomery(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = MRed(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryLvl multiplies p1 by p2 coefficient wise with a Montgomery
// modular reduction for the moduli q_0 up to q_level and returns the result on p3.
func (r *Ring) MulCoeffsMontgomeryLvl(level uint64, p1, p2, p3 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = MRed(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryAndAdd multiplies p1 by p2 coefficient wise with a
// Montgomery modular reduction and adds the result to p3.
func (r *Ring) MulCoeffsMontgomeryAndAdd(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = CRed(p3tmp[j]+MRed(p1tmp[j], p2tmp[j], qi, mredParams), qi)
		}
	}
}

// MulCoeffsMontgomeryAndAddLvl multiplies p1 by p2 coefficient wise with a Montgomery
// modular reduction for the moduli q_0 up to q_level and adds the result to p3.
func (r *Ring) MulCoeffsMontgomeryAndAddLvl(level uint64, p1, p2, p3 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = CRed(p3tmp[j]+MRed(p1tmp[j], p2tmp[j], qi, mredParams), qi)
		}
	}
}

// MulCoeffsMontgomeryAndAddNoMod multiplies p1 by p2 coefficient wise with a
// Montgomery modular reduction and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndAddNoMod(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] += MRed(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryAndAddNoModLvl multiplies p1 by p2 coefficient wise with a Montgomery modular
// reduction for the moduli q_0 up to q_level and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndAddNoModLvl(level uint64, p1, p2, p3 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] += MRed(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryConstantAndAddNoModLvl multiplies p1 by p2 coefficient wise with a constant time Montgomery
// modular reduction for the moduli q_0 up to q_level and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryConstantAndAddNoModLvl(level uint64, p1, p2, p3 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] += MRedConstant(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryAndSub multiplies p1 by p2 coefficient wise with
// a Montgomery modular reduction and subtracts the result to p3.
func (r *Ring) MulCoeffsMontgomeryAndSub(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = CRed(p3tmp[j]+(qi-MRed(p1tmp[j], p2tmp[j], qi, mredParams)), qi)
		}
	}
}

// MulCoeffsMontgomeryAndSubNoMod multiplies p1 by p2 coefficient wise with a Montgomery
// modular reduction and subtracts the result to p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndSubNoMod(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = p3tmp[j] + (qi - MRed(p1tmp[j], p2tmp[j], qi, mredParams))
		}
	}
}

// MulCoeffsConstant multiplies p1 by p2 coefficient wise with a constant
// time Barrett modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsConstant(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = BRedConstant(p1tmp[j], p2tmp[j], qi, bredParams)
		}
	}
}

// MulCoeffsMontgomeryConstant multiplies p1 by p2 coefficient wise with a
// constant time Montgomery modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsMontgomeryConstant(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p3tmp[j] = MRedConstant(p1tmp[j], p2tmp[j], qi, mredParams)
		}
	}
}

// MulPoly multiplies p1 by p2 and writes the result on p3.
func (r *Ring) MulPoly(p1, p2, p3 *Poly) {

	a := r.NewPoly()
	b := r.NewPoly()

	r.NTT(p1, a)
	r.NTT(p2, b)
	r.MulCoeffs(a, b, p3)
	r.InvNTT(p3, p3)
}

// MulPolyMontgomery multiplies p1 by p2 and writes the result on p3..
func (r *Ring) MulPolyMontgomery(p1, p2, p3 *Poly) {

	a := r.NewPoly()
	b := r.NewPoly()

	r.NTT(p1, a)
	r.NTT(p2, b)
	r.MulCoeffsMontgomery(a, b, p3)
	r.InvNTT(p3, p3)
}

// MulPolyNaive multiplies p1 by p2 with a naive convolution and writes the result on p3.
func (r *Ring) MulPolyNaive(p1, p2, p3 *Poly) {

	p1Copy := p1.CopyNew()
	p2Copy := p2.CopyNew()

	r.MForm(p1Copy, p1Copy)

	r.AND(p3, 0, p3)

	for x, qi := range r.Modulus {

		p1tmp, p2tmp, p3tmp := p1Copy.Coeffs[x], p2Copy.Coeffs[x], p3.Coeffs[x]

		mredParams := r.MredParams[x]

		for i := uint64(0); i < r.N; i++ {

			for j := uint64(0); j < i; j++ {
				p3tmp[j] = CRed(p3tmp[j]+(qi-MRed(p1tmp[i], p2tmp[r.N-i+j], qi, mredParams)), qi)
			}

			for j := uint64(i); j < r.N; j++ {
				p3tmp[j] = CRed(p3tmp[j]+MRed(p1tmp[i], p2tmp[j-i], qi, mredParams), qi)
			}
		}
	}
}

// MulPolyNaiveMontgomery multiplies p1 by p2 with a naive convolution and writes the result on p3.
func (r *Ring) MulPolyNaiveMontgomery(p1, p2, p3 *Poly) {

	p1Copy := p1.CopyNew()
	p2Copy := p2.CopyNew()

	r.AND(p3, 0, p3)

	for x, qi := range r.Modulus {

		p1tmp, p2tmp, p3tmp := p1Copy.Coeffs[x], p2Copy.Coeffs[x], p3.Coeffs[x]

		mredParams := r.MredParams[x]

		for i := uint64(0); i < r.N; i++ {

			for j := uint64(0); j < i; j++ {
				p3tmp[j] = CRed(p3tmp[j]+(qi-MRed(p1tmp[i], p2tmp[r.N-i+j], qi, mredParams)), qi)
			}

			for j := uint64(i); j < r.N; j++ {
				p3tmp[j] = CRed(p3tmp[j]+MRed(p1tmp[i], p2tmp[j-i], qi, mredParams), qi)
			}
		}
	}
}

// Exp raises p1 to p1^e and writes the result on p2.
// TODO : implement Montgomery ladder
func (r *Ring) Exp(p1 *Poly, e uint64, p2 *Poly) {

	r.NTT(p1, p1)

	tmp := r.NewPoly()
	r.Add(tmp, p1, tmp)

	for i := range r.Modulus {
		p2tmp := p2.Coeffs[i]
		for x := uint64(0); x < r.N; x++ {
			p2tmp[x] = 1
		}
	}

	for i := e; i > 0; i >>= 1 {
		if (i & 1) == 1 {
			r.MulCoeffs(p2, tmp, p2)
		}
		r.MulCoeffs(tmp, p1, tmp)
	}

	r.InvNTT(p2, p2)
	r.InvNTT(p1, p2)
}

// AddScalar adds to each coefficient of p1 a scalar and writes the result on p2.
func (r *Ring) AddScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	for i, Qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = CRed(p1tmp[j]+scalar, Qi)
		}
	}
}

// AddScalarBigint adds to each coefficient of p1 a big.Int scalar and writes the result on p2.
func (r *Ring) AddScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	tmp := new(big.Int)
	for i, Qi := range r.Modulus {
		scalarQi := tmp.Mod(scalar, NewUint(Qi)).Uint64()
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = CRed(p1tmp[j]+scalarQi, Qi)
		}
	}
}

// SubScalar subtracts to each coefficient of p1 a scalar and writes the result on p2.
func (r *Ring) SubScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	for i, Qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = CRed(p1tmp[j]+(Qi-scalar), Qi)
		}
	}
}

// SubScalarBigint subtracts to each coefficient of p1 a big.Int scalar and writes the result on p2.
func (r *Ring) SubScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	tmp := new(big.Int)
	for i, Qi := range r.Modulus {
		scalarQi := tmp.Mod(scalar, NewUint(Qi)).Uint64()
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = CRed(p1tmp[j]+(Qi-scalarQi), Qi)
		}
	}
}

// MulScalar multiplies each coefficient of p1 by a scalar and writes the result on p2.
func (r *Ring) MulScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	for i, Qi := range r.Modulus {
		scalarMont := MForm(BRedAdd(scalar, Qi, r.BredParams[i]), Qi, r.BredParams[i])
		mredParams := r.MredParams[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = MRed(p1tmp[j], scalarMont, Qi, mredParams)
		}
	}
}

// MulScalarLvl multiplies each coefficient of p1 by a scalar for the moduli q_0 up to q_level and writes the result on p2.
func (r *Ring) MulScalarLvl(level uint64, p1 *Poly, scalar uint64, p2 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		Qi := r.Modulus[i]
		scalarMont := MForm(BRedAdd(scalar, Qi, r.BredParams[i]), Qi, r.BredParams[i])
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = MRed(p1tmp[j], scalarMont, Qi, mredParams)
		}
	}
}

// MulScalarBigint multiplies each coefficientsof p1 by a big.Int scalar and writes the result on p2.
func (r *Ring) MulScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	scalarQi := new(big.Int)
	for i, Qi := range r.Modulus {
		scalarQi.Mod(scalar, NewUint(Qi))
		scalarMont := MForm(BRedAdd(scalarQi.Uint64(), Qi, r.BredParams[i]), Qi, r.BredParams[i])
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = MRed(p1tmp[j], scalarMont, Qi, mredParams)
		}
	}
}

// MulScalarBigintLvl multiplies each coefficientsof p1 by a big.Int scalar
//for the moduli q_0 up to q_level and writes the result on p2.
func (r *Ring) MulScalarBigintLvl(level uint64, p1 *Poly, scalar *big.Int, p2 *Poly) {
	scalarQi := new(big.Int)
	for i := uint64(0); i < level+1; i++ {
		Qi := r.Modulus[i]
		scalarQi.Mod(scalar, NewUint(Qi))
		scalarMont := MForm(BRedAdd(scalarQi.Uint64(), Qi, r.BredParams[i]), Qi, r.BredParams[i])
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = MRed(p1tmp[j], scalarMont, Qi, mredParams)
		}
	}
}

// Shift circulary shifts the coefficients of the polynomial p1 by n to the left and writes the result on p2.
func (r *Ring) Shift(p1 *Poly, n uint64, p2 *Poly) {
	mask := uint64((1 << r.N) - 1)
	for i := range r.Modulus {
		p2.Coeffs[i] = append(p1.Coeffs[i][(n&mask):], p1.Coeffs[i][:(n&mask)]...)
	}
}

// MForm switches p1 to the Montgomery domain and writes the result on p2.
func (r *Ring) MForm(p1, p2 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = MForm(p1tmp[j], qi, bredParams)
		}
	}
}

// MFormLvl switches p1 to the Montgomery domain for the moduli q_0 up to q_level and writes the result on p2.
func (r *Ring) MFormLvl(level uint64, p1, p2 *Poly) {
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		bredParams := r.BredParams[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = MForm(p1tmp[j], qi, bredParams)
		}
	}
}

// InvMForm switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
func (r *Ring) InvMForm(p1, p2 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = InvMForm(p1tmp[j], qi, mredParams)
		}
	}
}

// MulByPow2New multiplies p1 by 2^pow2 and returns the result a new polynomial p2.
func (r *Ring) MulByPow2New(p1 *Poly, pow2 uint64) (p2 *Poly) {
	p2 = r.NewPoly()
	r.MulByPow2(p1, pow2, p2)
	return
}

// MulByPow2 multiplies p1 by 2^pow2 and writes the result on p2.
func (r *Ring) MulByPow2(p1 *Poly, pow2 uint64, p2 *Poly) {
	r.MForm(p1, p2)
	for i, Qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = PowerOf2(p1tmp[j], pow2, Qi, mredParams)
		}
	}
}

// MulByPow2Lvl multiplies p1 by 2^pow2 for the moduli q_0 up to q_level and writes the result on p2.
func (r *Ring) MulByPow2Lvl(level uint64, p1 *Poly, pow2 uint64, p2 *Poly) {
	r.MFormLvl(level, p1, p2)
	for i := uint64(0); i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = PowerOf2(p1tmp[j], pow2, qi, mredParams)
		}
	}
}

// MultByMonomialNew multiplies p1 by x^monomialDeg and writes the result on a new polynomial p2.
func (r *Ring) MultByMonomialNew(p1 *Poly, monomialDeg uint64) (p2 *Poly) {
	p2 = r.NewPoly()
	r.MultByMonomial(p1, monomialDeg, p2)
	return
}

// MultByMonomial multiplies p1 by x^monomialDeg and writes the result on p2.
func (r *Ring) MultByMonomial(p1 *Poly, monomialDeg uint64, p2 *Poly) {

	shift := monomialDeg % (r.N << 1)

	if shift == 0 {

		for i := range r.Modulus {
			p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
			for j := uint64(0); j < r.N; j++ {
				p2tmp[j] = p1tmp[j]
			}
		}

	} else {

		tmpx := r.NewPoly()

		if shift < r.N {

			for i := range r.Modulus {
				p1tmp, tmpxT := p1.Coeffs[i], tmpx.Coeffs[i]
				for j := uint64(0); j < r.N; j++ {
					tmpxT[j] = p1tmp[j]
				}
			}

		} else {

			for i, qi := range r.Modulus {
				p1tmp, tmpxT := p1.Coeffs[i], tmpx.Coeffs[i]
				for j := uint64(0); j < r.N; j++ {
					tmpxT[j] = qi - p1tmp[j]
				}
			}
		}

		shift %= r.N

		for i, qi := range r.Modulus {
			p2tmp, tmpxT := p2.Coeffs[i], tmpx.Coeffs[i]
			for j := uint64(0); j < shift; j++ {
				p2tmp[j] = qi - tmpxT[r.N-shift+j]
			}
		}

		for i := range r.Modulus {
			p2tmp, tmpxT := p2.Coeffs[i], tmpx.Coeffs[i]
			for j := shift; j < r.N; j++ {
				p2tmp[j] = tmpxT[j-shift]

			}
		}
	}
}

// MulByVectorMontgomery multiplies p1 by a vector of uint64 coefficients and writes the result on p2.
func (r *Ring) MulByVectorMontgomery(p1 *Poly, vector []uint64, p2 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] = MRed(p1tmp[j], vector[j], qi, mredParams)
		}
	}
}

// MulByVectorMontgomeryAndAddNoMod multiplies p1 by a vector of uint64 coefficients and adds the result on p2 without modular reduction.
func (r *Ring) MulByVectorMontgomeryAndAddNoMod(p1 *Poly, vector []uint64, p2 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := uint64(0); j < r.N; j++ {
			p2tmp[j] += MRed(p1tmp[j], vector[j], qi, mredParams)
		}
	}
}

// BitReverse applies a bit reverse permutation on the coefficients of p1 and writes the result on p2.
// Can safely be used for inplace permutation.
func (r *Ring) BitReverse(p1, p2 *Poly) {
	bitLenOfN := uint64(bits.Len64(r.N) - 1)

	if p1 != p2 {
		for i := range r.Modulus {
			p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
			for j := uint64(0); j < r.N; j++ {
				p2tmp[utils.BitReverse64(j, bitLenOfN)] = p1tmp[j]
			}
		}
	} else { // In place in case p1 = p2
		for x := range r.Modulus {
			p2tmp := p2.Coeffs[x]
			for i := uint64(0); i < r.N; i++ {
				j := utils.BitReverse64(i, bitLenOfN)
				if i < j {
					p2tmp[i], p2tmp[j] = p2tmp[j], p2tmp[i]
				}
			}
		}
	}
}

// Rotate applies a Galoi Automorphism on p1 in NTT form,
// rotating the coefficients to the right by n and writes the result on p2.
// Requires the data to permuted in bitreversal order before applying NTT.
func (r *Ring) Rotate(p1 *Poly, n uint64, p2 *Poly) {

	var root, gal uint64

	n &= (1 << r.N) - 1

	for i, qi := range r.Modulus {

		mredParams := r.MredParams[i]

		root = MRed(r.PsiMont[i], r.PsiMont[i], qi, mredParams)

		root = modexpMontgomery(root, n, qi, mredParams, r.BredParams[i])

		gal = MForm(1, qi, r.BredParams[i])

		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]

		for j := uint64(1); j < r.N; j++ {

			gal = MRed(gal, root, qi, mredParams)

			p2tmp[j] = MRed(p1tmp[j], gal, qi, mredParams)

		}
	}
}
package ring

import (
	"github.com/ldsec/lattigo/utils"
	"math"
	"math/bits"
)

type TernarySampler struct {
	prng                    utils.PRNG
	context                 *Context
	matrixTernary           [][]uint64
	matrixTernaryMontgomery [][]uint64
}

// NewTernarySampler creates a new instance of TernarySampler.
// Accepts a PRNG and context and samples different kinds of ternary polynomials
func NewTernarySampler(prng utils.PRNG, context *Context) *TernarySampler {
	ternarySampler := new(TernarySampler)
	ternarySampler.context = context
	ternarySampler.prng = prng
	ternarySampler.InitialiseMatrix()
	return ternarySampler
}

func (ternarySampler *TernarySampler) InitialiseMatrix() {
	ternarySampler.matrixTernary = make([][]uint64, len(ternarySampler.context.Modulus))
	ternarySampler.matrixTernaryMontgomery = make([][]uint64, len(ternarySampler.context.Modulus))

	for i, Qi := range ternarySampler.context.Modulus {

		ternarySampler.matrixTernary[i] = make([]uint64, 3)
		ternarySampler.matrixTernary[i][0] = 0
		ternarySampler.matrixTernary[i][1] = 1
		ternarySampler.matrixTernary[i][2] = Qi - 1

		ternarySampler.matrixTernaryMontgomery[i] = make([]uint64, 3)
		ternarySampler.matrixTernaryMontgomery[i][0] = 0
		ternarySampler.matrixTernaryMontgomery[i][1] = MForm(1, Qi, ternarySampler.context.bredParams[i])
		ternarySampler.matrixTernaryMontgomery[i][2] = MForm(Qi-1, Qi, ternarySampler.context.bredParams[i])
	}
}

// SampleUniform samples a ternary polynomial with distribution [1/3, 1/3, 1/3].
func (ternarySampler *TernarySampler) SampleUniform(pol *Poly) {
	ternarySampler.sample(ternarySampler.matrixTernary, 1.0/3.0, pol)
}

// Sample samples a ternary polynomial with distribution [(1-p)/2, p, (1-p)/2].
func (ternarySampler *TernarySampler) Sample(pol *Poly, p float64) {
	ternarySampler.sample(ternarySampler.matrixTernary, p, pol)
}

// SampleMontgomery samples a ternary polynomial with distribution [(1-p)/2, p, (1-p)/2] in Montgomery form.
func (ternarySampler *TernarySampler) SampleMontgomery(pol *Poly, p float64) {
	ternarySampler.sample(ternarySampler.matrixTernaryMontgomery, p, pol)
}

// SampleNew samples a new ternary polynomial with distribution [(1-p)/2, p, (1-p)/2].
func (ternarySampler *TernarySampler) SampleNew(p float64) (pol *Poly) {
	pol = ternarySampler.context.NewPoly()
	ternarySampler.Sample(pol, p)
	return pol
}

// SampleMontgomeryNew samples a nes ternary polynomial with distribution [(1-p)/2, p, (1-p)/2] in Montgomery form.
func (ternarySampler *TernarySampler) SampleMontgomeryNew(p float64) (pol *Poly) {
	pol = ternarySampler.context.NewPoly()
	ternarySampler.SampleMontgomery(pol, p)
	return
}

// SampleNTTNew samples a new ternary polynomial with distribution [(1-p)/2, p, (1-p)/2] in NTT form.
func (ternarySampler *TernarySampler) SampleNTTNew(p float64) (pol *Poly) {
	pol = ternarySampler.context.NewPoly()
	ternarySampler.Sample(pol, p)
	ternarySampler.context.NTT(pol, pol)
	return
}

// SampleNTT samples a ternary polynomial with distribution [(1-p)/2, p, (1-p)/2] in NTT form.
func (ternarySampler *TernarySampler) SampleNTT(pol *Poly, p float64) {
	ternarySampler.Sample(pol, p)
	ternarySampler.context.NTT(pol, pol)
}

// SampleMontgomeryNTTNew samples a new ternary polynomial with distribution [(1-p)/2, p, (1-p)/2] in NTT and Montgomery form.
func (ternarySampler *TernarySampler) SampleMontgomeryNTTNew(p float64) (pol *Poly) {
	pol = ternarySampler.SampleMontgomeryNew(p)
	ternarySampler.context.NTT(pol, pol)
	return
}

// SampleMontgomeryNTT samples a ternary polynomial with distribution [(1-p)/2, p, (1-p)/2] in NTT and Montgomery form.
func (ternarySampler *TernarySampler) SampleMontgomeryNTT(pol *Poly, p float64) {
	ternarySampler.SampleMontgomery(pol, p)
	ternarySampler.context.NTT(pol, pol)
}

// SampleSparse samples a polynomial with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients.
func (ternarySampler *TernarySampler) SampleSparse(pol *Poly, hw uint64) {
	ternarySampler.sampleSparse(ternarySampler.matrixTernary, pol, hw)
}

// SampleSparseNew samples a new polynomial with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients.
func (ternarySampler *TernarySampler) SampleSparseNew(hw uint64) (pol *Poly) {
	pol = ternarySampler.context.NewPoly()
	ternarySampler.SampleSparse(pol, hw)
	return pol
}

// SampleSparseNTT samples a polynomial with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients in NTT form.
func (ternarySampler *TernarySampler) SampleSparseNTT(pol *Poly, hw uint64) {
	ternarySampler.SampleSparse(pol, hw)
	ternarySampler.context.NTT(pol, pol)
}

// SampleSparseNTTNew samples a new polynomial with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients in NTT form.
func (ternarySampler *TernarySampler) SampleSparseNTTNew(hw uint64) (pol *Poly) {
	pol = ternarySampler.context.NewPoly()
	ternarySampler.sampleSparse(ternarySampler.matrixTernaryMontgomery, pol, hw)
	ternarySampler.context.NTT(pol, pol)
	return pol
}

// SampleSparseMontgomery samples a polynomial with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients in Montgomery form.
func (ternarySampler *TernarySampler) SampleSparseMontgomery(pol *Poly, hw uint64) {
	ternarySampler.sampleSparse(ternarySampler.matrixTernaryMontgomery, pol, hw)
}

// SampleSparseMontgomeryNew samples a new polynomial with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients in Montgomery form.
func (ternarySampler *TernarySampler) SampleSparseMontgomeryNew(hw uint64) (pol *Poly) {
	pol = ternarySampler.context.NewPoly()
	ternarySampler.sampleSparse(ternarySampler.matrixTernaryMontgomery, pol, hw)
	return pol
}

// SampleSparseMontgomeryNTTNew samples a new polynomial with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients in NTT and Montgomery form.
func (ternarySampler *TernarySampler) SampleSparseMontgomeryNTTNew(hw uint64) (pol *Poly) {
	pol = ternarySampler.SampleSparseMontgomeryNew(hw)
	ternarySampler.context.NTT(pol, pol)
	return pol
}

// SampleSparseMontgomeryNTT samples a polynomial with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients in NTT and Montgomery form.
func (ternarySampler *TernarySampler) SampleSparseMontgomeryNTT(pol *Poly, hw uint64) {
	ternarySampler.SampleSparseMontgomery(pol, hw)
	ternarySampler.context.NTT(pol, pol)
}

func computeMatrixTernary(p float64) (M [][]uint8) {
	var g float64
	var x uint64

	precision := uint64(56)

	M = make([][]uint8, 2)

	g = p
	g *= math.Exp2(float64(precision))
	x = uint64(g)

	M[0] = make([]uint8, precision-1)
	for j := uint64(0); j < precision-1; j++ {
		M[0][j] = uint8((x >> (precision - j - 1)) & 1)
	}

	g = 1 - p
	g *= math.Exp2(float64(precision))
	x = uint64(g)

	M[1] = make([]uint8, precision-1)
	for j := uint64(0); j < precision-1; j++ {
		M[1][j] = uint8((x >> (precision - j - 1)) & 1)
	}

	return M
}

func (ternarySampler *TernarySampler) sample(samplerMatrix [][]uint64, p float64, pol *Poly) {

	if p == 0 {
		panic("cannot sample -> p = 0")
	}

	var coeff uint64
	var sign uint64
	var index uint64

	if p == 0.5 {

		randomBytesCoeffs := make([]byte, ternarySampler.context.N>>3)
		randomBytesSign := make([]byte, ternarySampler.context.N>>3)

		ternarySampler.prng.Clock(randomBytesCoeffs)

		ternarySampler.prng.Clock(randomBytesSign)

		for i := uint64(0); i < ternarySampler.context.N; i++ {
			coeff = uint64(uint8(randomBytesCoeffs[i>>3])>>(i&7)) & 1
			sign = uint64(uint8(randomBytesSign[i>>3])>>(i&7)) & 1

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j := range ternarySampler.context.Modulus {
				pol.Coeffs[j][i] = samplerMatrix[j][index] //(coeff & (sign^1)) | (qi - 1) * (sign & coeff)
			}
		}

	} else {

		matrix := computeMatrixTernary(p)

		randomBytes := make([]byte, ternarySampler.context.N)

		pointer := uint8(0)
		bytePointer := uint64(0)

		ternarySampler.prng.Clock(randomBytes)

		for i := uint64(0); i < ternarySampler.context.N; i++ {

			coeff, sign, randomBytes, pointer, bytePointer = kysampling(ternarySampler.prng, matrix, randomBytes, pointer, bytePointer, ternarySampler.context.N)

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j := range ternarySampler.context.Modulus {
				pol.Coeffs[j][i] = samplerMatrix[j][index] //(coeff & (sign^1)) | (qi - 1) * (sign & coeff)
			}
		}
	}
}

func (ternarySampler *TernarySampler) sampleSparse(samplerMatrix [][]uint64, pol *Poly, hw uint64) {

	if hw > ternarySampler.context.N {
		hw = ternarySampler.context.N
	}

	var mask, j uint64
	var coeff uint8

	index := make([]uint64, ternarySampler.context.N)
	for i := uint64(0); i < ternarySampler.context.N; i++ {
		index[i] = i
	}

	randomBytes := make([]byte, (uint64(math.Ceil(float64(hw) / 8.0)))) // We sample ceil(hw/8) bytes
	pointer := uint8(0)

	ternarySampler.prng.Clock(randomBytes)

	for i := uint64(0); i < hw; i++ {
		mask = (1 << uint64(bits.Len64(ternarySampler.context.N-i))) - 1 // rejection sampling of a random variable between [0, len(index)]

		j = randInt32(ternarySampler.prng, mask)
		for j >= ternarySampler.context.N-i {
			j = randInt32(ternarySampler.prng, mask)
		}

		coeff = (uint8(randomBytes[0]) >> (i & 7)) & 1 // random binary digit [0, 1] from the random bytes
		for i := range ternarySampler.context.Modulus {
			pol.Coeffs[i][index[j]] = samplerMatrix[i][coeff]
		}

		// Removes the element in position j of the slice (order not preserved)
		index[j] = index[len(index)-1]
		index = index[:len(index)-1]

		pointer++

		if pointer == 8 {
			randomBytes = randomBytes[1:]
			pointer = 0
		}
	}
}

// kysampling use the binary expension and random bytes matrix to sample a discret gaussian value and its sign.
func kysampling(prng utils.PRNG, M [][]uint8, randomBytes []byte, pointer uint8, bytePointer uint64, byteLength uint64) (uint64, uint64, []byte, uint8, uint64) {

	var sign uint8

	d := 0
	col := 0
	colLen := len(M)

	for {

		// Uses one random byte per cycle and cycle through the randombytes
		for i := pointer; i < 8; i++ {

			d = (d << 1) + 1 - int((uint8(randomBytes[bytePointer])>>i)&1)

			// There is small probability that it will get out of the bound, then
			// rerun until it gets a proper output
			if d > colLen-1 {
				return kysampling(prng, M, randomBytes, i, bytePointer, byteLength)
			}

			for row := colLen - 1; row >= 0; row-- {

				d -= int(M[row][col])

				if d == -1 {

					// Sign
					if i == 7 {
						pointer = 0
						// If the last bit of the array was read, samples a new one
						bytePointer++

						if bytePointer >= byteLength {
							bytePointer = 0
							prng.Clock(randomBytes)
						}

						sign = uint8(randomBytes[bytePointer]) & 1

					} else {
						pointer = i
						// Else the sign is the next bit of the byte
						sign = uint8(randomBytes[bytePointer]>>(i+1)) & 1
					}

					return uint64(row), uint64(sign), randomBytes, pointer + 1, bytePointer
				}
			}

			col++
		}

		// Resets the bit pointer and discards the used byte
		pointer = 0
		// If the last bit of the array was read, samples a new one
		bytePointer++

		if bytePointer >= byteLength {
			bytePointer = 0
			prng.Clock(randomBytes)
		}

	}
}
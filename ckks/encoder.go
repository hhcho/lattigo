package ckks

import (
	"fmt"
	"github.com/hhcho/mpc-core"
	"github.com/ldsec/lattigo/ring"
	"math"
	"math/big"
	"math/rand"
	"github.com/hhcho/frand"

)

// Encoder is an interface implenting the encoding algorithms.
type Encoder interface {
	Encode(plaintext *Plaintext, values []complex128, slots uint64)
	EncodeNew(values []complex128, slots uint64) (plaintext *Plaintext)
	Decode(plaintext *Plaintext, slots uint64) (res []complex128)

	EncodeRVecNew(values mpc_core.RVec, slots uint64, fracBits int) (plaintext *Plaintext)
	DecodeRVec(plaintext *Plaintext, slots uint64, fracBits int) (res mpc_core.RVec)

	FFTTest(values []complex128, N uint64, prg *frand.RNG)
}

// encoder is a struct storing the necessary parameters to encode a slice of complex number on a Plaintext.
type encoder struct {
	params       *Parameters
	ckksContext  *Context
	values       []complex128
	valuesfloat  []float64
	bigintCoeffs []*big.Int
	qHalf        *big.Int
	polypool     *ring.Poly
	m            uint64
	roots        []complex128
	rootsBig     []bigComplex
	rotGroup     []uint64
}

type bigComplex struct {
	real *big.Float
	imag *big.Float
}

func (encoder *encoder) FFTTest(values []complex128, slots uint64, prg *frand.RNG) {
	fmt.Println("FFTTest")
	//cos test
	for i := 0.0; i < 2*3.1415; i += 0.2 {
		f, _ := cosBig(big.NewFloat(i)).Float64()
		fmt.Println("cos(",i,"): ", f, math.Cos(i), f-math.Cos(i))
		f, _ = sinBig(big.NewFloat(i)).Float64()
		fmt.Println("sin(",i,"): ", f, math.Sin(i),  f-math.Sin(i))
	}
	return


	fmt.Println("Modulus:", encoder.ckksContext.contextQ.Modulus)
	fmt.Println("input:", values[:5])
	encoder.invfft(values, slots)
	fmt.Println("invfft:", values[:5])

	fmt.Println("maxSlots:", encoder.ckksContext.maxSlots, "slots:", slots)
	fmt.Println("scale:", encoder.params.Scale)
	fmt.Println("log2(scale):", math.Log2(encoder.params.Scale))

	vec := make([]uint32, slots)
	vecfloat := make([]float64, slots)
	for i := range vec {
		vec[i] = rand.Uint32()
		vecfloat[i] = float64(vec[i]) / encoder.params.Scale
	}

	for i, jdx, idx := uint64(0), encoder.ckksContext.maxSlots, uint64(0); i < slots; i, jdx, idx = i+1, jdx+1, idx+1 {
		encoder.valuesfloat[idx] = real(encoder.values[i])
		encoder.valuesfloat[jdx] = imag(encoder.values[i])
	}

	fmt.Println("len(valuesfloat):", len(encoder.valuesfloat))
	fmt.Println("valuesfloat:", encoder.valuesfloat[:5])

	rtype := mpc_core.LElem128Zero
	mask := make(mpc_core.RVec, len(encoder.valuesfloat))
	for i := range mask {
		mask[i] = rtype.RandBits(prg, rtype.ModBitLength() - 4)

		if i < 5 {
			fmt.Println(i, ":", mask[i].(mpc_core.LElem128).ToBigInt().String())
		}
	}

	pt := NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.Scale)
	pt2 := NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.Scale)

	scaleUpVecExact(vecfloat, pt.scale, encoder.ckksContext.contextQ.Modulus[:pt.Level()+1], pt.value.Coeffs)

	moduli := encoder.ckksContext.contextQ.Modulus[:pt.Level()+1]
	var xInt *big.Int
	tmp := new(big.Int)
	for i := range mask {
		//xInt = new(big.Int)
		//xInt.SetUint64(uint64(vec[i]))
		xInt = mask[i].(mpc_core.LElem128).ToBigInt()

		for j := range moduli {
			tmp.Mod(xInt, ring.NewUint(moduli[j]))
			pt2.value.Coeffs[j][i] = tmp.Uint64()
		}
	}

	for j := range pt.value.Coeffs {
		fmt.Println(j, ":", pt.value.Coeffs[j][:5])
		fmt.Println(j, ":", pt2.value.Coeffs[j][:5])
	}

	//encoder.ckksContext.contextQ.NTTLvl(pt.Level(), pt.value, pt.value)

	//encoder.fft(values, slots)
	//fmt.Println("fft:", values[:5])
}

var piBase, _ = new(big.Float).SetPrec(200).SetString("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127372458700660631558817488152092096282925409171536436789259036001133053054882046652138414695194151160943305727036575959195309218611738193261179310511854807446237996274956735188575272489122793818301194912")
var piTable = make(map[uint]*big.Float)
func piBig(prec uint) *big.Float {
	if prec > piBase.Prec() {
		//fmt.Println("Warning: requested precision for Pi ", prec, " exceeds constant ", piBase.Prec())
		prec = piBase.Prec()
	}

	v, flag := piTable[prec]
	if flag {
		return v
	}

	fmt.Println("Pi", prec)

	v = new(big.Float)
	v.Copy(piBase)
	v.SetPrec(prec)
	piTable[prec] = v

	return v
}
func sinBig(x *big.Float) *big.Float {
	tmp := new(big.Float)
	tmp.Mul(piBig(x.Prec()), big.NewFloat(0.5))
	tmp.Sub(x, tmp)
	return cosBig(tmp) // six(x) = cos(x - pi/2)
}
func cosBig(x *big.Float) *big.Float {
	x = new(big.Float).Abs(x)
	tmp := new(big.Float)
	pi := piBig(x.Prec())
	flag2pi := x.Cmp(tmp.Mul(pi, big.NewFloat(2)))
	flagpi := x.Cmp(pi)
	flagpi2 := x.Cmp(tmp.Mul(pi, big.NewFloat(0.5)))
	if flag2pi > 0 {
		panic("cosBig input outside of [-2*pi, 2*pi]")
	} else if flag2pi == 0 {
		return big.NewFloat(1)
	} else if flagpi >= 0 {
		return tmp.Neg(cosBig(tmp.Sub(x, pi)))
	} else if flagpi2 == 0 {
		return big.NewFloat(0)
	} else if flagpi2 > 0 {
		return tmp.Neg(cosBig(tmp.Sub(x, pi)))
	}

	// x is now in [0, pi/2]

	// Choose the number of steps k of the Cordic algorithm
	// k=250 gives about 150 decimals, k=1000 about 600
	steps := 250 // TODO: choose based on the precision of x

	// xˆ2/2ˆ(2k)
	t := new(big.Float)
	t.SetMantExp(big.NewFloat(1), -2*steps)
	s := new(big.Float)
	tmp.Mul(x, x)
	s.Mul(tmp, t)

	four := big.NewFloat(4)
	for i := 1; i <= steps; i++ {
		tmp.Sub(four, s)
		s.Mul(s, tmp)
	}

	out := new(big.Float)
	tmp.Quo(s, big.NewFloat(2))
	out.Sub(big.NewFloat(1), tmp)
	return out
}

// NewEncoder creates a new Encoder that is used to encode a slice of complex values of size at most N/2 (the number of slots) on a Plaintext.
func NewEncoder(params *Parameters) Encoder {

	if !params.isValid {
		panic("cannot newEncoder: parameters are invalid (check if the generation was done properly)")
	}

	m := uint64(2 << params.LogN)

	rotGroup := make([]uint64, m>>1)
	fivePows := uint64(1)
	for i := uint64(0); i < m>>2; i++ {
		rotGroup[i] = fivePows
		fivePows *= GaloisGen
		fivePows &= (m - 1)
	}

	var angle float64
	roots := make([]complex128, m+1)
	rootsBig := make([]bigComplex, m+1)
	// TODO replace hardcoded precision
	anglePart := new(big.Float).Mul(piBig(128), big.NewFloat(2))
	anglePart.Quo(anglePart, big.NewFloat(float64(m)))

	for i := uint64(0); i < m; i++ {
		angle = 2 * 3.141592653589793 * float64(i) / float64(m)
		roots[i] = complex(math.Cos(angle), math.Sin(angle))
		angleBig := new(big.Float)
		angleBig.Mul(anglePart, big.NewFloat(float64(i)))
		rootsBig[i] = bigComplex{cosBig(angleBig), sinBig(angleBig)}
	}
	roots[m] = roots[0]
	rootsBig[m] = rootsBig[0]

	ckksContext := newContext(params)

	return &encoder{
		params:       params.Copy(),
		ckksContext:  ckksContext,
		values:       make([]complex128, m>>2),
		valuesfloat:  make([]float64, m>>1),
		bigintCoeffs: make([]*big.Int, m>>1),
		qHalf:        ring.NewUint(0),
		polypool:     ckksContext.contextQ.NewPoly(),
		m:            m,
		rotGroup:     rotGroup,
		roots:        roots,
		rootsBig:	  rootsBig,
	}
}

func (encoder *encoder) EncodeNew(values []complex128, slots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.Scale)
	encoder.Encode(plaintext, values, slots)
	return
}

// Encode takes a slice of complex128 values of size at most N/2 (the number of slots) and encodes it in the receiver Plaintext.
func (encoder *encoder) Encode(plaintext *Plaintext, values []complex128, slots uint64) {

	if uint64(len(values)) > encoder.ckksContext.maxSlots || uint64(len(values)) > slots {
		panic("cannot Encode: too many values for the given number of slots")
	}

	if slots == 0 && slots&(slots-1) == 0 {
		panic("cannot Encode: slots must be a power of two between 1 and N/2")
	}

	if uint64(len(values)) != slots {
		panic("cannot Encode: number of values must be equal to slots")
	}

	for i := uint64(0); i < slots; i++ {
		encoder.values[i] = values[i]
	}

	encoder.invfft(encoder.values, slots)

	gap := encoder.ckksContext.maxSlots / slots

	for i, jdx, idx := uint64(0), encoder.ckksContext.maxSlots, uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		encoder.valuesfloat[idx] = real(encoder.values[i])
		encoder.valuesfloat[jdx] = imag(encoder.values[i])
	}

	scaleUpVecExact(encoder.valuesfloat, plaintext.scale, encoder.ckksContext.contextQ.Modulus[:plaintext.Level()+1], plaintext.value.Coeffs)

	encoder.ckksContext.contextQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)

	for i := uint64(0); i < encoder.ckksContext.maxSlots; i++ {
		encoder.values[i] = 0
	}

	for i := uint64(0); i < encoder.ckksContext.n; i++ {
		encoder.valuesfloat[i] = 0
	}
}

func (encoder *encoder) EncodeRVecNew(values mpc_core.RVec, slots uint64, fracBits int) (plaintext *Plaintext) {
	if uint64(len(values)) > encoder.ckksContext.maxSlots || uint64(len(values)) > slots {
		panic("cannot EncodeRVecNew: too many values for the given number of slots")
	}

	if slots == 0 {
		return NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.Scale)
	}

	rtype := mpc_core.LElem128Zero

	if slots&(slots-1) != 0 { // slots not a power of two
		closestPow := uint64(math.Pow(2, math.Ceil(math.Log2(float64(slots)))))
		newValues := mpc_core.InitRVec(rtype.Zero(), int(closestPow))
		for i := range values {
			newValues[i] = values[i]
		}
		values = newValues
		slots = closestPow
	}

	if uint64(len(values)) != slots {
		panic("cannot EncodeRVecNew: number of values must be equal to slots")
	}

	if values.Type().TypeID() != rtype.TypeID() {
		panic("cannot EncodeRVecNew: only LElem128 supported")
	}

	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.Scale)

	slice := make([]bigComplex, len(encoder.values))
	zeroFloat := big.NewFloat(0)
	for i := range slice {
		if uint64(i) < slots {
			slice[i] = bigComplex{values[i].(mpc_core.LElem128).ToSignedBigFloat(fracBits), big.NewFloat(0)}
		} else {
			slice[i] = bigComplex{zeroFloat, zeroFloat}
		}
	}

	encoder.invfftBig(slice, slots)

	gap := encoder.ckksContext.maxSlots / slots

	floatSlice := make([]*big.Float, len(encoder.valuesfloat))
	for i := range floatSlice {
		floatSlice[i] = zeroFloat
	}

	for i, jdx, idx := uint64(0), encoder.ckksContext.maxSlots, uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		floatSlice[idx] = slice[i].real
		floatSlice[jdx] = slice[i].imag
	}

	moduli := encoder.ckksContext.contextQ.Modulus[:plaintext.Level()+1]
	xInt := new(big.Int)
	xFlo := new(big.Float)
	scaleBig := big.NewFloat(plaintext.scale)
	tmp := new(big.Int)
	for i := range floatSlice {
		xFlo.Mul(floatSlice[i], scaleBig)
		xInt = arithRound(xFlo)

		for j := range moduli {
			tmp.Mod(xInt, ring.NewUint(moduli[j]))
			plaintext.value.Coeffs[j][i] = tmp.Uint64()
		}
	}

	encoder.ckksContext.contextQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)
	return
}


func (encoder *encoder) DecodeRVec(plaintext *Plaintext, slots uint64, fracBits int) (res mpc_core.RVec) {
	rtype := mpc_core.LElem128Zero

	encoder.ckksContext.contextQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	encoder.ckksContext.contextQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

	Q := encoder.ckksContext.bigintChain[plaintext.Level()]

	maxSlots := encoder.ckksContext.maxSlots

	encoder.qHalf.Set(Q)
	encoder.qHalf.Rsh(encoder.qHalf, 1)

	gap := encoder.ckksContext.maxSlots / slots

	values := make([]bigComplex, len(encoder.values))

	var sign int

	for i, idx := uint64(0), uint64(0); i < slots; i, idx = i+1, idx+gap {

		// Centers the value around the current modulus
		encoder.bigintCoeffs[idx].Mod(encoder.bigintCoeffs[idx], Q)
		sign = encoder.bigintCoeffs[idx].Cmp(encoder.qHalf)
		if sign == 1 || sign == 0 {
			encoder.bigintCoeffs[idx].Sub(encoder.bigintCoeffs[idx], Q)
		}

		// Centers the value around the current modulus
		encoder.bigintCoeffs[idx+maxSlots].Mod(encoder.bigintCoeffs[idx+maxSlots], Q)
		sign = encoder.bigintCoeffs[idx+maxSlots].Cmp(encoder.qHalf)
		if sign == 1 || sign == 0 {
			encoder.bigintCoeffs[idx+maxSlots].Sub(encoder.bigintCoeffs[idx+maxSlots], Q)
		}

		values[i] = bigComplex{scaleDownBig(encoder.bigintCoeffs[idx], plaintext.scale), scaleDownBig(encoder.bigintCoeffs[idx+maxSlots], plaintext.scale)}
	}

	encoder.fftBig(values, slots)

	res = mpc_core.InitRVec(rtype.Zero(), int(slots))

	tmpFloat := new(big.Float)
	scale := big.NewFloat(0)
	scale.SetInt(new(big.Int).Lsh(big.NewInt(1), uint(fracBits)))
	for i := uint64(0); i < slots; i++ {
		v := values[i].real
		if v.Sign() < 0 {
			tmpFloat.Mul(tmpFloat.Neg(v), scale)
			res[i] = res[i].Sub(rtype.FromBigInt(arithRound(tmpFloat)))
		} else {
			tmpFloat.Mul(v, scale)
			res[i] = rtype.FromBigInt(arithRound(tmpFloat))
		}
	}

	return res
}

// Decode decodes the Plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *encoder) Decode(plaintext *Plaintext, slots uint64) (res []complex128) {

	encoder.ckksContext.contextQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	encoder.ckksContext.contextQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

	Q := encoder.ckksContext.bigintChain[plaintext.Level()]

	maxSlots := encoder.ckksContext.maxSlots

	encoder.qHalf.Set(Q)
	encoder.qHalf.Rsh(encoder.qHalf, 1)

	gap := encoder.ckksContext.maxSlots / slots

	var sign int

	for i, idx := uint64(0), uint64(0); i < slots; i, idx = i+1, idx+gap {

		// Centers the value around the current modulus
		encoder.bigintCoeffs[idx].Mod(encoder.bigintCoeffs[idx], Q)
		sign = encoder.bigintCoeffs[idx].Cmp(encoder.qHalf)
		if sign == 1 || sign == 0 {
			encoder.bigintCoeffs[idx].Sub(encoder.bigintCoeffs[idx], Q)
		}

		// Centers the value around the current modulus
		encoder.bigintCoeffs[idx+maxSlots].Mod(encoder.bigintCoeffs[idx+maxSlots], Q)
		sign = encoder.bigintCoeffs[idx+maxSlots].Cmp(encoder.qHalf)
		if sign == 1 || sign == 0 {
			encoder.bigintCoeffs[idx+maxSlots].Sub(encoder.bigintCoeffs[idx+maxSlots], Q)
		}

		encoder.values[i] = complex(scaleDown(encoder.bigintCoeffs[idx], plaintext.scale), scaleDown(encoder.bigintCoeffs[idx+maxSlots], plaintext.scale))
	}

	encoder.fft(encoder.values, slots)

	res = make([]complex128, slots)

	for i := range res {
		res[i] = encoder.values[i]

	}

	for i := uint64(0); i < encoder.ckksContext.maxSlots; i++ {
		encoder.values[i] = 0
	}

	return
}

func (encoder *encoder) invfftlazy(values []complex128, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v complex128

	for len := N; len >= 1; len >>= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (lenq - (encoder.rotGroup[j] % lenq)) * gap
				u = values[i+j] + values[i+j+lenh]
				v = values[i+j] - values[i+j+lenh]
				v *= encoder.roots[idx]
				values[i+j] = u
				values[i+j+lenh] = v

			}
		}
	}

	sliceBitReverseInPlaceComplex128(values, N)
}

func (encoder *encoder) invfft(values []complex128, N uint64) {

	encoder.invfftlazy(values, N)

	for i := uint64(0); i < N; i++ {
		values[i] /= complex(float64(N), 0)
	}
}

func (encoder *encoder) fft(values []complex128, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v complex128

	sliceBitReverseInPlaceComplex128(values, N)

	for len := uint64(2); len <= N; len <<= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (encoder.rotGroup[j] % lenq) * gap
				u = values[i+j]
				v = values[i+j+lenh]
				v *= encoder.roots[idx]
				values[i+j] = u + v
				values[i+j+lenh] = u - v
			}
		}
	}
}


func (encoder *encoder) invfftlazyBig(values []bigComplex, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v bigComplex

	for len := N; len >= 1; len >>= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (lenq - (encoder.rotGroup[j] % lenq)) * gap
				//u = values[i+j] + values[i+j+lenh]
				u = values[i+j].copy()
				u.add(values[i+j+lenh])

				//v = values[i+j] - values[i+j+lenh]
				v = values[i+j].copy()
				v.sub(values[i+j+lenh])

				//v *= encoder.roots[idx]
				v.mul(encoder.rootsBig[idx])

				values[i+j] = u
				values[i+j+lenh] = v

			}
		}
	}

	sliceBitReverseInPlaceBigComplex(values, N)
}

func (encoder *encoder) invfftBig(values []bigComplex, N uint64) {

	encoder.invfftlazyBig(values, N)

	for i := uint64(0); i < N; i++ {
		//values[i] /= complex(float64(N), 0)
		values[i].real.Quo(values[i].real, big.NewFloat(float64(N)))
		values[i].imag.Quo(values[i].imag, big.NewFloat(float64(N)))
	}
}

func (encoder *encoder) fftBig(values []bigComplex, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v bigComplex

	sliceBitReverseInPlaceBigComplex(values, N)

	for len := uint64(2); len <= N; len <<= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (encoder.rotGroup[j] % lenq) * gap
				u = values[i+j]
				v = values[i+j+lenh].copy()

				//v *= encoder.roots[idx]
				v.mul(encoder.rootsBig[idx])

				//values[i+j] = u + v
				values[i+j] = u.copy()
				values[i+j].add(v)

				//values[i+j+lenh] = u - v
				values[i+j+lenh] = u.copy()
				values[i+j+lenh].sub(v)
			}
		}
	}
}

func (bc bigComplex) copy() (out bigComplex) {
	out.imag = new(big.Float)
	out.real = new(big.Float)
	out.imag.Copy(bc.imag)
	out.real.Copy(bc.real)
	return out
}

func (bc *bigComplex) fromComplex128(v complex128) {
	bc.imag = big.NewFloat(imag(v))
	bc.real = big.NewFloat(real(v))
}

func (bc *bigComplex) add(v bigComplex) {
	bc.real.Add(bc.real, v.real)
	bc.imag.Add(bc.imag, v.imag)
}

func (bc *bigComplex) sub(v bigComplex) {
	bc.real.Sub(bc.real, v.real)
	bc.imag.Sub(bc.imag, v.imag)
}

func (bc *bigComplex) mul(v bigComplex) {
	t1 := new(big.Float)
	t2 := new(big.Float)
	outr := new(big.Float)
	outi := new(big.Float)

	t1.Mul(bc.real, v.real)
	t2.Mul(bc.imag, v.imag)
	outr.Sub(t1, t2)

	t1.Mul(bc.real, v.imag)
	t2.Mul(bc.imag, v.real)
	outi.Add(t1, t2)

	bc.real = outr
	bc.imag = outi
}

func sliceBitReverseInPlaceBigComplex(slice []bigComplex, N uint64) {

	var bit, j uint64

	for i := uint64(1); i < N; i++ {

		bit = N >> 1

		for j >= bit {
			j -= bit
			bit >>= 1
		}

		j += bit

		if i < j {
			slice[i], slice[j] = slice[j], slice[i]
		}
	}
}

func scaleDownBig(coeff *big.Int, n float64) (x *big.Float) {

	x = new(big.Float).SetInt(coeff)
	x.Quo(x, big.NewFloat(n))

	return
}

var halfFloat = big.NewFloat(0.5)
func arithRound(a *big.Float) *big.Int {
	var i *big.Int
	if a.Signbit() {
		i, _ = new(big.Float).Sub(a, halfFloat).Int(nil)
	} else {
		i, _ = new(big.Float).Add(a, halfFloat).Int(nil)
	}
	return i
}
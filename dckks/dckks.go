package dckks

import (
	"github.com/ldsec/lattigo/ckks"
	"github.com/ldsec/lattigo/ring"
	"math"
)

type Context struct {
	params *ckks.Parameters

	n uint64

	gaussianSampler *ring.KYSampler

	ContextQ  *ring.Context
	ContextP  *ring.Context
	ContextQP *ring.Context

	Alpha uint64
	Beta  uint64
}

func NewContext(params *ckks.Parameters) (context *Context) {

	if !params.IsValid() {
		panic("cannot NewContext : params not valid (check if they where generated properly)")
	}

	context = new(Context)

	context.params = params.Copy()

	n := uint64(1 << params.LogN)

	context.n = n

	context.Alpha = uint64(len(params.Pi))
	context.Beta = uint64(math.Ceil(float64(len(params.Qi)) / float64(context.Alpha)))

	var err error
	if context.ContextQ, err = ring.NewContextWithParams(n, params.Qi); err != nil {
		panic(err)
	}

	if context.ContextP, err = ring.NewContextWithParams(n, params.Pi); err != nil {
		panic(err)
	}

	if context.ContextQP, err = ring.NewContextWithParams(n, append(params.Qi, params.Pi...)); err != nil {
		panic(err)
	}

	context.gaussianSampler = context.ContextQP.NewKYSampler(params.Sigma, int(params.Sigma*6))

	return
}

// NewCRPGenerator creates a CRPGenerator
func NewCRPGenerator(params *ckks.Parameters, key []byte) *ring.CRPGenerator {
	ctx := NewContext(params)
	return ring.NewCRPGenerator(key, ctx.ContextQP)
}

package quadratures

import "math"

type Quadrature struct {
	nodes   []float64
	weights []float64
	f       func(float64) float64
}

func (q *Quadrature) Integrate(a, b, h float64) float64 {
	result := 0.0
	numIntervals := int((b - a) / h)

	for interval := 0; interval < numIntervals; interval++ {
		subIntervalStart := a + float64(interval)*h

		for node := 0; node < len(q.nodes); node++ {
			x := (h/2.0)*(1+q.nodes[node]) + subIntervalStart
			result += q.weights[node] * q.f(x)
		}
	}

	return result * (h / 2.0)
}

func NewGauss4(f func(x float64) float64) *Quadrature {
	return &Quadrature{
		nodes: []float64{ // Узлы: ±sqrt((3 ± 2*sqrt(6/5)) / 7)
			-math.Sqrt((3 + 2*math.Sqrt(6.0/5.0)) / 7.0),
			-math.Sqrt((3 - 2*math.Sqrt(6.0/5.0)) / 7.0),
			math.Sqrt((3 - 2*math.Sqrt(6.0/5.0)) / 7.0),
			math.Sqrt((3 + 2*math.Sqrt(6.0/5.0)) / 7.0),
		},

		weights: []float64{ // Веса: (18 ± sqrt(30)) / 36
			(18 - math.Sqrt(30.0)) / 36.0,
			(18 + math.Sqrt(30.0)) / 36.0,
			(18 + math.Sqrt(30.0)) / 36.0,
			(18 - math.Sqrt(30.0)) / 36.0,
		},

		f: f,
	}
}

func NewSimpson(f func(x float64) float64) *Quadrature {
	return &Quadrature{
		nodes: []float64{
			-1,
			0,
			1,
		},

		weights: []float64{
			1.0 / 3.0,
			4.0 / 3.0,
			1.0 / 3.0,
		},

		f: f,
	}
}

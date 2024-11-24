package main

import (
	"cm5/quadratures"
	"fmt"
	"math"
)

func f(x float64) float64 {
	return math.Tan(x) + x
}

func testMethod(method *quadratures.Quadrature, a, b, starI float64) {
	h := math.Abs(b - a)

	for i := 0; i <= 2; i++ {
		I := method.Integrate(a, b, h)

		halfI := method.Integrate(a, b, h/2)
		k := math.Log2(math.Abs(1 + (halfI-I)/(starI-halfI)))

		rungeError := (halfI - I) / (math.Pow(2, k) - 1)

		fmt.Printf("%-30s%e\n", "h:", h)
		fmt.Printf("%-30s%e\n", "I*-Ih:", starI-I)
		fmt.Printf("%-30s%e\n", "(I*-Ih) / (I*-Ih2):", (starI-I)/(starI-halfI))
		fmt.Printf("%-30s%e\n", "(Ih/2-Ih) / (2^k - 1):", rungeError)
		fmt.Printf("%-30s%e\n", "Ir:", halfI+rungeError)
		fmt.Printf("%-30s%e\n", "I* - Ir:", starI-halfI-rungeError)
		fmt.Printf("%-30s%e\n", "k:", k)

		fmt.Printf("___________________________________________________________________________________\n")

		h /= 2
	}
}

func main() {
	a := 0.0
	b := math.Pi / 4
	starI := 0.654998727814

	method := quadratures.NewSimpson(f)
	fmt.Println("Simpson method testing:")
	testMethod(method, a, b, starI)

	fmt.Println("\n\n")

	method = quadratures.NewGauss4(f)
	fmt.Println("Gauss-4 method testing:")
	testMethod(method, a, b, starI)

}

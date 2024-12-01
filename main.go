package main

import (
	"encoding/json"
	"fmt"
	"os"
	"sync"
)

func main() {
	var wg sync.WaitGroup
	// wg.Add(1)
	// go problem1(&wg)
	// wg.Add(1)
	// go problem2(&wg)
	// wg.Add(1)
	// go problem3(&wg)
	wg.Add(1)
	go problem7(&wg, 0.1)
	wg.Add(1)
	go problem7(&wg, 1)
	wg.Add(1)
	go problem7(&wg, 10)
	wg.Wait()
}

func problem1(w *sync.WaitGroup) {
	defer w.Done()
	ratios := []float64{0.001, 0.003, 0.01, 0.03, 0.1, 0.3}
	mus := make([]*Mu, len(ratios))
	msds := make([][]float64, len(ratios))
	times := 1000
	particles := 1000
	phi := 0.05
	sigma := 1.0
	l := L(sigma, particles, phi)
	var m sync.Mutex
	var wg sync.WaitGroup
	for i, ratio := range ratios {
		wg.Add(1)
		go func() {
			defer wg.Done()
			defer fmt.Printf("Finished ratio = %f\n", ratio)
			mu := NewMu(sigma, particles, phi, ratio, l, l, 0)
			msd := []float64{}
			for i := 0; i < times; i++ {
				mu.Evolve()
				msd = append(msd, mu.MSD())
			}
			m.Lock()
			msds[i] = msd
			mus[i] = mu
			m.Unlock()
		}()
	}
	wg.Wait()

	output := map[string]interface{}{
		"times": times,
		"msds":  msds,
		"mus":   mus,
	}

	SavAsJson(output, "output1")
}

func problem2(w *sync.WaitGroup) {
	defer w.Done()
	ratio := 0.1
	msds := [][]float64{}
	times := 1000
	particles := 1000
	phis := []float64{0.05, 0.2, 0.5}
	sigma := 1.0
	mus := []*Mu{}
	var wg sync.WaitGroup
	for _, phi := range phis {
		l := L(sigma, particles, phi)
		wg.Add(1)
		go func() {
			defer wg.Done()
			defer fmt.Printf("Finished phi = %f\n", phi)
			mu := NewMuWithTrianglarLattice(sigma, particles, phi, ratio, int(l), int(l), 0)
			msd := []float64{}
			for i := 0; i < times; i++ {
				mu.EvolveWithTriangle()
				msd = append(msd, mu.MSD())
			}
			msds = append(msds, msd)
			mus = append(mus, mu)
		}()
	}
	wg.Wait()

	output := map[string]interface{}{
		"times": times,
		"msds":  msds,
		"mus":   mus,
	}
	SavAsJson(output, "problem2/output")
}
func problem3(w *sync.WaitGroup) {
	defer w.Done()
	ratio := 0.1
	times := 100000
	particles := 1000
	gravities := []float64{0, 0.01, 0.1, 1, 10}
	// gravities := []float64{0.1}
	sigma := 1.0
	mus := []*Mu{}
	phi := 0.05
	l := L(sigma, particles, phi)
	t := 1.0
	var wg sync.WaitGroup
	for _, gravity := range gravities {
		wg.Add(1)
		go func() {
			defer wg.Done()
			defer fmt.Printf("Finished gravity = %f\n", gravity)
			mu := NewMuWithWall(sigma, particles, phi, ratio, l, l*10, gravity, t)
			for i := 0; i < times; i++ {
				mu.EvolveWithWall()
			}
			mus = append(mus, mu)
		}()
	}
	wg.Wait()

	output := map[string]interface{}{
		"times": times,
		"mus":   mus,
	}
	SavAsJson(output, "problem3/output_test")
}

func problem7(w *sync.WaitGroup, t float64) {
	defer w.Done()
	ratio := 0.1
	times := 100
	particles := 1000
	gravities := []float64{0, 0.01, 0.1, 1, 10}
	// tempertures := []float64{0.1, 1, 10}
	// tempertures := []float64{1}
	// gravities := []float64{0.1}
	// t := 1.0
	sigma := 1.0
	mus := []*Mu{}
	phi := 0.05
	l := L(sigma, particles, phi)
	var wg sync.WaitGroup
	for _, gravity := range gravities {
		wg.Add(1)
		go func() {
			defer wg.Done()
			defer fmt.Printf("Finished temp: %v, gravity = %f\n", t, gravity)
			mu := NewMuWithWall(sigma, particles, phi, ratio, l, l*10, gravity, t)
			for i := 0; i < times; i++ {
				mu.EvolveWithWall()
			}
			mus = append(mus, mu)
		}()
	}
	wg.Wait()

	output := map[string]interface{}{
		"times": times,
		"mus":   mus,
	}
	SavAsJson(output, fmt.Sprintf("problem7/output_%v", t))
}

func SavAsJson(output interface{}, filename string) {
	file, err := os.Create("./values/" + filename + ".json")
	if err != nil {
		fmt.Println("Error creating file:", err)
		return
	}
	defer file.Close()

	encoder := json.NewEncoder(file)
	encoder.SetIndent("", "  ")
	if err := encoder.Encode(output); err != nil {
		fmt.Println("Error encoding JSON:", err)
		return
	}
	fmt.Printf("a file %v saved\n", filename)
}

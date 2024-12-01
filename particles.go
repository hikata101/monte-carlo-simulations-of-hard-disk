package main

import (
	"fmt"
	"math"
	"math/rand"
)

func NewParticle(x, y float64, sigma float64) *Particle {
	return &Particle{
		InitialX: x,
		InitialY: y,
		X:        x,
		Y:        y,
		Sigma:    sigma,
	}
}

// Particle
type Particle struct {
	InitialX float64
	InitialY float64
	X        float64
	Y        float64
	Sigma    float64
}

func (p *Particle) Configurations() {
	fmt.Println(p.InitialX, p.InitialY, p.X, p.Y, p.Sigma)
}

// Calculate square displacement
func (p *Particle) SquareDisplacement() float64 {
	return (p.X-p.InitialX)*(p.X-p.InitialX) + (p.Y-p.InitialY)*(p.Y-p.InitialY)
}

// configuration
func L(sigma float64, n int, phi float64) float64 {
	return math.Sqrt(math.Pi * sigma * sigma * float64(n) / 4 / phi)
}

func NewMu(sigma float64, n int, phi float64, deltaOverSigma float64, lx, ly float64, gravity float64) *Mu {
	mu := &Mu{
		Particles: []*Particle{},
		Sigma:     sigma,
		Lx:        lx,
		Ly:        ly,
		Phi:       phi,
		Gravity:   gravity,
		Delta:     sigma * deltaOverSigma,
		Ratio:     deltaOverSigma,
		Kb:        1,
		T:         1,
		M:         1,
	}
	for i := 0; i < n; i++ {
		signx := 1.0
		if rand.Float64() < 0.5 {
			signx = -1.0
		}
		x := rand.Float64() * lx / 2 * signx
		signy := 1.0
		if rand.Float64() < 0.5 {
			signy = -1.0
		}
		y := rand.Float64() * ly / 2 * signy
		if mu.Accept(i, x, y) {
			mu.Particles = append(mu.Particles, NewParticle(x, y, sigma))
		} else {
			i--
		}
	}
	fmt.Printf("initialization of mu has done. sigma: %v, Lx: %v, Ly: %v \n", mu.Sigma, mu.Lx, mu.Ly)
	return mu
}

func NewMuWithTrianglarLattice(sigma float64, n int, phi float64, deltaOverSigma float64, lx, ly int, gravity float64) *Mu {
	mu := &Mu{
		Particles: []*Particle{},
		Sigma:     sigma,
		Lx:        float64(lx),
		Ly:        float64(ly),
		Phi:       phi,
		Gravity:   gravity,
		Delta:     sigma * deltaOverSigma,
		Ratio:     deltaOverSigma,
		Kb:        1,
		T:         1,
		M:         1,
	}
	for i := 0; i < n; i++ {
		xIndex := rand.Intn(lx)
		x := float64(xIndex) * sigma * 1.5
		y := float64(rand.Intn(ly)+(xIndex%2)) * sigma * 1.5 * math.Sqrt(3) / 2
		if mu.Accept(i, x, y) {
			mu.Particles = append(mu.Particles, NewParticle(x, y, sigma))
		} else {
			i--
		}
	}
	fmt.Printf("initialization of mu has done. sigma: %v, Lx: %v, Ly: %v \n", mu.Sigma, mu.Lx, mu.Ly)
	return mu
}

func NewMuWithTrianglarLatticeAndWall(sigma float64, n int, phi float64, deltaOverSigma float64, lx, ly int, gravity float64) *Mu {
	mu := &Mu{
		Particles: []*Particle{},
		Sigma:     sigma,
		Lx:        float64(lx),
		Ly:        float64(ly),
		Phi:       phi,
		Gravity:   gravity,
		Delta:     sigma * deltaOverSigma,
		Ratio:     deltaOverSigma,
		Kb:        1,
		T:         1,
		M:         1,
	}
	for i := 0; i < n; i++ {
		xIndex := rand.Intn(lx)
		x := float64(xIndex) * sigma * 1.5
		y := float64(rand.Intn(ly)+(xIndex%2)) * sigma * 1.5 * math.Sqrt(3) / 2
		if !WallPositive(float64(lx), float64(ly), x, y) {
			i--
			continue
		}
		if mu.Accept(i, x, y) {
			mu.Particles = append(mu.Particles, NewParticle(x, y, sigma))
		} else {
			i--
		}
	}
	fmt.Printf("initialization of mu has done. sigma: %v, Lx: %v, Ly: %v \n", mu.Sigma, mu.Lx, mu.Ly)
	return mu
}

func NewMuWithWall(sigma float64, n int, phi float64, deltaOverSigma float64, lx, ly float64, gravity, t float64) *Mu {
	mu := &Mu{
		Particles: []*Particle{},
		Sigma:     sigma,
		Lx:        lx,
		Ly:        ly,
		Phi:       phi,
		Gravity:   gravity,
		Delta:     sigma * deltaOverSigma,
		Ratio:     deltaOverSigma,
		Kb:        1,
		T:         t,
		M:         1,
	}
	for i := 0; i < n; i++ {
		signx := 1.0
		if rand.Float64() < 0.5 {
			signx = -1.0
		}
		x := rand.Float64() * lx / 2 * signx
		signy := 1.0
		if rand.Float64() < 0.5 {
			signy = -1.0
		}
		y := rand.Float64() * ly / 2 * signy
		if !Wall(lx, ly, x, y, sigma) {
			i--
			continue
		}
		if !mu.Collide(x, y) {
			mu.Particles = append(mu.Particles, NewParticle(x, y, sigma))
		} else {
			i--
		}
	}
	fmt.Printf("initialization of mu has done. sigma: %v, Lx: %v, Ly: %v, g: %v, t: %v \n", mu.Sigma, mu.Lx, mu.Ly, mu.Gravity, mu.T)
	return mu
}

// Mu is a collection of particles
type Mu struct {
	Particles []*Particle
	Sigma     float64 // particle radius
	Delta     float64 // delta = sigma * deltaOverSigma
	Ratio     float64 // ratio of delta
	Lx        float64 // boundary length x
	Ly        float64 // boundary length y
	Phi       float64 // area fraction
	Gravity   float64 // gravity
	Accepted  int
	Denied    int
	Kb        float64 // Boltzmann constant = 1
	T         float64 // temperature
	M         float64 // mass
}

func (m *Mu) Collide(x, y float64) bool {
	for _, particle := range m.Particles {
		distance := math.Sqrt((x-particle.X)*(x-particle.X) + (y-particle.Y)*(y-particle.Y))
		if distance <= m.Sigma {
			return true
		}
	}
	return false
}

// Mean square displacement
func (mu *Mu) MSD() float64 {
	msd := 0.0
	for _, particle := range mu.Particles {
		msd += particle.SquareDisplacement()
	}
	return msd / float64(len(mu.Particles))
}

// Evolve all particles.
func (mu *Mu) EvolveWithTriangle() {
	for index, particle := range mu.Particles {
		xIndex := rand.Intn(int(mu.Lx))
		xEvolved := float64(xIndex) * mu.Sigma * 1.5
		yEvolved := float64(rand.Intn(int(mu.Ly))+(xIndex%2)) * mu.Sigma * 1.5 * math.Sqrt(3) / 2
		if mu.Accept(index, xEvolved, yEvolved) {
			mu.Accepted++
			particle.X = xEvolved
			particle.Y = yEvolved
		} else {
			mu.Denied++
		}
	}
}

// Evolve all particles with periodic boundary condition.
func (mu *Mu) Evolve() {
	for index, particle := range mu.Particles {
		xEvolved := EvolveValue(particle.X, mu.Delta, mu.Lx)
		yEvolved := EvolveValue(particle.Y, mu.Delta, mu.Ly)
		if mu.Accept(index, xEvolved, yEvolved) {
			mu.Accepted++
			particle.X = xEvolved
			particle.Y = yEvolved
		} else {
			mu.Denied++
		}
	}
}

// Evolve all particles with wall.
func (mu *Mu) EvolveWithWall() {
	for index, particle := range mu.Particles {
		xEvolved := EvolveWithWall(particle.X, mu.Delta, mu.Lx)
		yEvolved := EvolveWithWall(particle.Y, mu.Delta, mu.Ly)
		if mu.Accept(index, xEvolved, yEvolved) {
			mu.Accepted++
			particle.X = xEvolved
			particle.Y = yEvolved
		} else {
			mu.Denied++
		}
	}
}

func WallPositive(lx, ly, x, y float64) bool {
	if x > lx || x < 0 {
		return false
	}
	if y > ly || x < 0 {
		return false
	}
	return true
}

func Wall(lx, ly, x, y, sigma float64) bool {
	if math.Abs(lx/2-math.Abs(x)) < sigma {
		return false
	}
	if math.Abs(ly/2-math.Abs(y)) < sigma {
		return false
	}
	return true
}

func (mu *Mu) EvolveStepWithTriangle() {
	// a := mu.Sigma // 格子間隔と粒子の直径を同じに設定

	// 三角格子の6つの方向ベクトルを定義
	// directions := [][2]float64{
	// 	{1 * a, 0},                        // 右
	// 	{-1 * a, 0},                       // 左
	// 	{0.5 * a, math.Sqrt(3) / 2 * a},   // 右上
	// 	{-0.5 * a, math.Sqrt(3) / 2 * a},  // 左上
	// 	{0.5 * a, -math.Sqrt(3) / 2 * a},  // 右下
	// 	{-0.5 * a, -math.Sqrt(3) / 2 * a}, // 左下
	// }

	// 全粒子を時間発展
	for i := range mu.Particles {
		dx := mu.Sigma * 1.5
		dy := mu.Sigma * 1.5 * math.Sqrt(3) / 2
		// ランダムな方向を選択
		// dir := directions[rand.Intn(len(directions))]

		// 新しい位置の計算
		// xNew := mu.Particles[i].X + dir[0]
		// yNew := mu.Particles[i].Y + dir[1]
		xNew := mu.Particles[i].X + dx*Sign()
		yNew := mu.Particles[i].Y + dy*Sign()
		// 吸収なし、反射なし
		if !WallPositive(xNew, yNew, mu.Lx, mu.Ly) {
			continue
		}
		// // エネルギー変化を計算
		// dE := 0.0
		// for j, other := range mu.Particles {
		// 	if i == j {
		// 		continue
		// 	}
		// 	// エネルギー差を計算（現位置と新位置のエネルギー差）
		// 	dE += U(mu.Sigma, other.X, other.Y, xNew, yNew, mu.M, mu.Gravity) - U(mu.Sigma, other.X, other.Y, mu.Particles[i].X, mu.Particles[i].Y, mu.M, mu.Gravity)
		// }

		// Metropolis法による受け入れ判定
		if mu.Accept(i, xNew, yNew) {
			mu.Accepted++
			mu.Particles[i].X = xNew
			mu.Particles[i].Y = yNew
		} else {
			mu.Denied++
		}
	}
}

func Sign() float64 {
	if rand.Float64() < 0.5 {
		return -1
	}
	return 1
}

// accept the trial move or not as the system's energy
func (mu *Mu) Accept(index int, xNew, yNew float64) bool {
	for i, other := range mu.Particles {
		if i == index {
			// Skip the particle itself
			continue
		}
		acceptance := 0.0
		if index > len(mu.Particles)-1 {
			energy := U(mu.Sigma, other.X, other.Y, xNew, yNew, mu.M, mu.Gravity)
			acceptance = Acceptance(math.Exp(-1 / mu.Kb / mu.T * energy))
			if acceptance == 0 {
				// Not Accept the trial move
				return false
			} else {
				continue
			}
		}
		// Calculate the energy difference
		dE := U(mu.Sigma, other.X, other.Y, xNew, yNew, mu.M, mu.Gravity) - U(mu.Sigma, other.X, other.Y, mu.Particles[index].X, mu.Particles[index].Y, mu.M, mu.Gravity)
		acceptanceProb := math.Exp(-dE / (mu.Kb * mu.T))
		if rand.Float64() >= acceptanceProb {
			return false
		}
	}
	return true
}

func (mu *Mu) Positions() ([]float64, []float64) {
	x := make([]float64, len(mu.Particles))
	y := make([]float64, len(mu.Particles))
	for i, particle := range mu.Particles {
		x[i] = particle.X
		y[i] = particle.Y
	}
	return x, y
}

func EvolveValue(value, delta float64, l float64) float64 {
	tmp := value + delta*(rand.Float64()-0.5)
	// Periodic boundary condition
	if tmp > l/2 {
		tmp -= l
	} else if tmp < -l/2 {
		tmp += l
	}
	return tmp
}

func EvolveWithWall(value, delta float64, l float64) float64 {
	tmp := value + delta*(rand.Float64()-0.5)
	if tmp > l/2 {
		return value
	} else if tmp < -l/2 {
		return value
	}
	return tmp
}

func PoseWallWithTriangle(x, y, lx, ly float64) bool {
	if math.Abs(x) > lx/2 || math.Abs(y) > ly/2 {
		return false
	}
	return true
}

// If r >= sigma, return 0. Otherwise, return infinity(Collision).
func U(sigma, x1, y1, x2, y2 float64, m, g float64) float64 {
	r := math.Sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1))
	if r >= sigma {
		return m * g * (y2 - y1)
		// return 0
	}
	return math.Inf(1)
}

func Ugravity(sigma, x1, y1, x2, y2, g, m, t, kb float64) float64 {
	r := math.Sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1))
	if r >= sigma {
		return 0
	}
	return 0.5 * m * g * r * r
}

// Dicide whether to accept the trial move or not.
// if probability is smaller than 1, return 0 which means accepted.
func Acceptance(probability float64) float64 {
	return math.Min(1, probability)
}

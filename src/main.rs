#[macro_use]
extern crate crunchy;
use std::ops::{Add, Sub, Mul, AddAssign, SubAssign};
use std::f64::consts::PI;

const SOLAR_MASS: f64 = 4.0 * PI * PI;
const YEAR: f64 = 365.24;

#[derive(Clone, Copy)]
struct Vec3(pub f64, pub f64, pub f64, f64);

impl Vec3 {
    fn new() -> Self {
        Vec3(0.0, 0.0, 0.0, 0.0)
    }

    fn norm(&self) -> f64 {
        self.squared_norm().sqrt()
    }

    fn squared_norm(&self) -> f64 {
        self.0 * self.0 + self.1 * self.1 + self.2 * self.2
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Vec3(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2, 0.0)
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Vec3(self.0 - rhs.0, self.1 - rhs.1, self.2 - rhs.2, 0.0)
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        *self = Vec3(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2, 0.0);
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = Vec3(self.0 - rhs.0, self.1 - rhs.1, self.2 - rhs.2, 0.0);
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Vec3(self.0 * rhs, self.1 * rhs, self.2 * rhs, 0.0)
    }
}

struct Planet {
    position: Vec3,
    velocity: Vec3,
    mass: f64,
}

impl Planet {
    fn sun() -> Self {
        Planet {
            position: Vec3(0.0, 0.0, 0.0, 0.0),
            velocity: Vec3(0.0, 0.0, 0.0, 0.0),
            mass: SOLAR_MASS,
        }
    }

    fn jupiter() -> Self {
        Planet {
            position: Vec3(
                4.84143144246472090e+00,
                -1.16032004402742839e+00,
                -1.03622044471123109e-01,
                0.0,
            ),
            velocity: Vec3(
                1.66007664274403694e-03 * YEAR,
                7.69901118419740425e-03 * YEAR,
                -6.90460016972063023e-05 * YEAR,
                0.0,
            ),
            mass: 9.54791938424326609e-04 * SOLAR_MASS,
        }
    }

    fn saturn() -> Self {
        Planet {
            position: Vec3(
                8.34336671824457987e+00,
                4.12479856412430479e+00,
                -4.03523417114321381e-01,
                0.0,
            ),
            velocity: Vec3(
                -2.76742510726862411e-03 * YEAR,
                4.99852801234917238e-03 * YEAR,
                2.30417297573763929e-05 * YEAR,
                0.0,
            ),
            mass: 2.85885980666130812e-04 * SOLAR_MASS,
        }
    }

    fn uranus() -> Self {
        Planet {
            position: Vec3(
                1.28943695621391310e+01,
                -1.51111514016986312e+01,
                -2.23307578892655734e-01,
                0.0,
            ),
            velocity: Vec3(
                2.96460137564761618e-03 * YEAR,
                2.37847173959480950e-03 * YEAR,
                -2.96589568540237556e-05 * YEAR,
                0.0,
            ),
            mass: 4.36624404335156298e-05 * SOLAR_MASS,
        }
    }

    fn neptune() -> Self {
        Planet {
            position: Vec3(
                1.53796971148509165e+01,
                -2.59193146099879641e+01,
                1.79258772950371181e-01,
                0.0,
            ),
            velocity: Vec3(
                2.68067772490389322e-03 * YEAR,
                1.62824170038242295e-03 * YEAR,
                -9.51592254519715870e-05 * YEAR,
                0.0,
            ),
            mass: 5.15138902046611451e-05 * SOLAR_MASS,
        }
    }
}

fn offset_momentum(velocity: &mut [Vec3], mass: &[f64]) {
    let total_momentum: Vec3 = velocity.iter().zip(mass.iter()).fold(Vec3::new(), |p,
     (&vel, &ma)| {
        p + (vel * ma)
    });
    velocity[0] = total_momentum * (-1.0 / mass[0]);
}

fn energy(position: &[Vec3], velocity: &[Vec3], mass: &[f64]) -> f64 {
    let mut energy = 0.0;
    for i in 0..5 {
        energy += velocity[i].squared_norm() * mass[i] / 2.0 -
            mass[i] *
                position
                    .iter()
                    .skip(i + 1)
                    .zip(mass.iter().skip(i + 1))
                    .map(|(pos, mass)| mass / (*pos - position[i]).norm())
                    .sum::<f64>()
    }
    energy
}

fn advance(position: &mut [Vec3], velocity: &mut [Vec3], mass: &[f64], dt: f64) {
    let mut r = [Vec3::new(); 10];
    let mut mag = [0.0; 10];

    let mut m = 0;
    unroll! {
        for i in 0..5 {
            for j in i + 1..5 {
                r[m] = position[i] - position[j];
                m += 1;
            }
        }
    }

    unroll!{
        for m in 0..10 {
            let d2 = r[m].squared_norm();
            let mut distance = 1.0 / d2.sqrt();
            distance = distance * (1.5 - 0.5 * d2 * distance * distance);
            mag[m] = dt * distance.powi(3);
        }
    }

    m = 0;
    unroll! {
        for i in 0..5 {
            for j in i + 1..5 {
                velocity[i] -= r[m] * mass[j] * mag[m];
                velocity[j] += r[m] * mass[i] * mag[m];
                m += 1;
            }
        }
    }

    unroll!{
        for m in 0..5 {
            position[m] += velocity[m] * dt;
        }
    }
}

fn main() {
    let n = std::env::args_os()
        .nth(1)
        .and_then(|s| s.into_string().ok())
        .and_then(|n| n.parse().ok())
        .unwrap_or(1000);

    let mut position = [
        Planet::sun().position,
        Planet::jupiter().position,
        Planet::saturn().position,
        Planet::uranus().position,
        Planet::neptune().position,
    ];
    let mut velocity = [
        Planet::sun().velocity,
        Planet::jupiter().velocity,
        Planet::saturn().velocity,
        Planet::uranus().velocity,
        Planet::neptune().velocity,
    ];
    let mass = [
        Planet::sun().mass,
        Planet::jupiter().mass,
        Planet::saturn().mass,
        Planet::uranus().mass,
        Planet::neptune().mass,
    ];

    offset_momentum(&mut velocity, &mass);
    println!("{:.9}", energy(&position, &velocity, &mass));
    for _ in 0..n {
        advance(&mut position, &mut velocity, &mass, 0.01);
    }
    println!("{:.9}", energy(&position, &velocity, &mass));
}
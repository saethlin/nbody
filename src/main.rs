#[macro_use]
extern crate crunchy;

use std::ops::{Add, Sub, Mul, AddAssign, SubAssign};
use std::f64::consts::PI;

const SOLAR_MASS: f64 = 4.0 * PI * PI;
const YEAR: f64 = 365.24;
const N_BODIES: usize = 5;
const N: usize = 10;

#[derive(Clone, Copy)]
struct Vec3(pub f64, pub f64, pub f64, pub f64);

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
    let total_momentum: Vec3 = velocity.iter().zip(mass.iter()).fold(
        Vec3::new(),
        |v, (vel, ma)| {
            v + *vel * *ma
        },
    );
    velocity[0] = total_momentum * (1.0 / mass[0]);
}

fn energy(position: &[Vec3], velocity: &[Vec3], mass: &[f64]) -> f64 {
    let mut energy = 0.0;
    for i in 0..N_BODIES {
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
    let i = 0;
    for j in 1..5 {
        r[m] = position[i] - position[j];
        m += 1;
    }

    let i = 1;
    for j in 2..5 {
        r[m] = position[i] - position[j];
        m += 1;
    }

    let i = 2;
    for j in 3..5 {
        r[m] = position[i] - position[j];
        m += 1;
    }

    let i = 3;
    for j in 4..5 {
        r[m] = position[i] - position[j];
        m += 1;
    }

    let i = 4;
    for j in 5..5 {
        r[m] = position[i] - position[j];
        m += 1;
    }
    unroll!{
    for m in 0..10 {
        let d2 = r[m].0.powi(2) + r[m].1.powi(2) + r[m].2.powi(2);
        let mut distance = 1.0 / d2.sqrt();
        distance = distance * (1.5 - 0.5 * d2 * distance * distance);
        mag[m] = dt * distance.powi(3);        
    }
    }
    /*
    let d2 = r[0].0.powi(2) + r[0].1.powi(2) + r[0].2.powi(2);
    let mut distance = 1.0 / d2.sqrt();
    distance = distance * (1.5 - 0.5 * d2 * distance * distance);
    mag[0] = dt * distance.powi(3);

    let d2 = r[1].0.powi(2) + r[1].1.powi(2) + r[1].2.powi(2);
    let mut distance = 1.0 / d2.sqrt();
    distance = distance * (1.5 - 0.5 * d2 * distance * distance);
    mag[1] = dt * distance.powi(3);

    let d2 = r[2].0.powi(2) + r[2].1.powi(2) + r[2].2.powi(2);
    let mut distance = 1.0 / d2.sqrt();
    distance = distance * (1.5 - 0.5 * d2 * distance * distance);
    mag[2] = dt * distance.powi(3);

    let d2 = r[3].0.powi(2) + r[3].1.powi(2) + r[3].2.powi(2);
    let mut distance = 1.0 / d2.sqrt();
    distance = distance * (1.5 - 0.5 * d2 * distance * distance);
    mag[3] = dt * distance.powi(3);

    let d2 = r[4].0.powi(2) + r[4].1.powi(2) + r[4].2.powi(2);
    let mut distance = 1.0 / d2.sqrt();
    distance = distance * (1.5 - 0.5 * d2 * distance * distance);
    mag[4] = dt * distance.powi(3);

    let d2 = r[5].0.powi(2) + r[5].1.powi(2) + r[5].2.powi(2);
    let mut distance = 1.0 / d2.sqrt();
    distance = distance * (1.5 - 0.5 * d2 * distance * distance);
    mag[5] = dt * distance.powi(3);

    let d2 = r[6].0.powi(2) + r[6].1.powi(2) + r[6].2.powi(2);
    let mut distance = 1.0 / d2.sqrt();
    distance = distance * (1.5 - 0.5 * d2 * distance * distance);
    mag[6] = dt * distance.powi(3);

    let d2 = r[7].0.powi(2) + r[7].1.powi(2) + r[7].2.powi(2);
    let mut distance = 1.0 / d2.sqrt();
    distance = distance * (1.5 - 0.5 * d2 * distance * distance);
    mag[7] = dt * distance.powi(3);

    let d2 = r[8].0.powi(2) + r[8].1.powi(2) + r[8].2.powi(2);
    let mut distance = 1.0 / d2.sqrt();
    distance = distance * (1.5 - 0.5 * d2 * distance * distance);
    mag[8] = dt * distance.powi(3);

    let d2 = r[9].0.powi(2) + r[9].1.powi(2) + r[9].2.powi(2);
    let mut distance = 1.0 / d2.sqrt();
    distance = distance * (1.5 - 0.5 * d2 * distance * distance);
    mag[9] = dt * distance.powi(3);
    */

    m = 0;
    unroll!{
    for i in 0..5 {
        for j in 1 + i..5 {
            let tmp = r[m] * mag[m];
            velocity[i] -= tmp * mass[j];
            velocity[j] += tmp * mass[i];
            m += 1;
        }
    }
    }

    /*
    m = 0;
    let i = 0;
    for j in 1..5 {
        velocity[i] -= r[m] * mass[j] * mag[m];
        velocity[j] += r[m] * mass[i] * mag[m];
        m += 1;
    }
    let i = 1;
    for j in 2..5 {
        velocity[i] -= r[m] * mass[j] * mag[m];
        velocity[j] += r[m] * mass[i] * mag[m];
        m += 1;
    }
    let i = 2;
    for j in 3..5 {
        velocity[i] -= r[m] * mass[j] * mag[m];
        velocity[j] += r[m] * mass[i] * mag[m];
        m += 1;
    }
    let i = 3;
    for j in 4..5 {
        velocity[i] -= r[m] * mass[j] * mag[m];
        velocity[j] += r[m] * mass[i] * mag[m];
        m += 1;
    }
    let i = 4;
    for j in 5..5 {
        velocity[i] -= r[m] * mass[j] * mag[m];
        velocity[j] += r[m] * mass[i] * mag[m];
        m += 1;
    }
    */

    position[0] += velocity[0] * dt;
    position[1] += velocity[1] * dt;
    position[2] += velocity[2] * dt;
    position[3] += velocity[3] * dt;
    position[4] += velocity[4] * dt;
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
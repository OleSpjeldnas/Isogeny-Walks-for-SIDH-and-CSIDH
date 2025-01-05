use super::F;

// Solves a linear system of equations
pub fn solve_linear_system(coeffs: Vec<Vec<F>>, constants: Vec<F>) -> Vec<F> {
    let size = coeffs.len();
    let (lower, upper) = lu_decomposition(&coeffs, size);

    let mut solutions = constants.clone();
    for i in 0..size {
        let mut sum = F::from(0);
        for j in 0..i {
            sum += lower[i][j] * solutions[j];
        }
        solutions[i] = (constants[i] - sum) / lower[i][i];
    }

    for i in (0..size).rev() {
        let mut sum = F::from(0);
        for j in (i + 1)..size {
            sum += upper[i][j] * solutions[j];
        }
        solutions[i] = (solutions[i] - sum) / upper[i][i];
    }

    solutions
}
// Computes the LU decomposition of a matrix
pub fn lu_decomposition(mat: &Vec<Vec<F>>, n: usize) -> (Vec<Vec<F>>, Vec<Vec<F>>) {
    let mut lower = vec![vec![F::from(0); n]; n];
    let mut upper = vec![vec![F::from(0); n]; n];

    for i in 0..n {
        for k in i..n {
            let mut sum = F::from(0);
            for j in 0..i {
                sum += lower[i][j] * upper[j][k];
            }
            upper[i][k] = mat[i][k] - sum;
        }

        for k in i..n {
            if i == k {
                lower[i][i] = F::from(1);
            } else {
                let mut sum = F::from(0);
                for j in 0..i {
                    sum += lower[k][j] * upper[j][i];
                }
                lower[k][i] = (mat[k][i] - sum) / upper[i][i];
            }
        }
    }

    (lower, upper)
}

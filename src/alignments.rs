// Return edit distance
pub fn edit_dist(x1: &[u8], x2: &[u8]) -> u32 {
    let mut dp = vec![vec![0; x2.len() + 1]; x1.len() + 1];
    for (i, row) in dp.iter_mut().enumerate() {
        row[0] = i as u32;
    }
    for j in 0..=x2.len() {
        dp[0][j] = j as u32;
    }
    for (i, x1_b) in x1.iter().enumerate() {
        for (j, x2_b) in x2.iter().enumerate() {
            let m = (x1_b != x2_b) as u32;
            dp[i + 1][j + 1] = (dp[i][j + 1] + 1).min(dp[i + 1][j] + 1).min(dp[i][j] + m);
        }
    }
    dp[x1.len()][x2.len()]
}

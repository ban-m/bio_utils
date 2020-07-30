// Return edit distance
pub fn edit_dist(x1: &[u8], x2: &[u8]) -> u32 {
    let mut dp = vec![vec![0; x2.len() + 1]; x1.len() + 1];
    for i in 0..=x1.len() {
        dp[i][0] = i as u32;
    }
    for j in 0..=x2.len() {
        dp[0][j] = j as u32;
    }
    for i in 0..x1.len() {
        for j in 0..x2.len() {
            let m = if x1[i] == x2[j] { 0 } else { 1 };
            dp[i + 1][j + 1] = (dp[i][j + 1] + 1).min(dp[i + 1][j] + 1).min(dp[i][j] + m);
        }
    }
    dp[x1.len()][x2.len()]
}

/// Iterates over a sequence slice with a sliding window,
/// calculates the average gc-content inside each window,
/// and returns an iterator over (window_start_index, average_gc_content)
pub fn get_gc_content(
    sequence: &[u8],
    window_size: usize,
) -> impl Iterator<Item = (usize, f32)> + '_ {
    let step = window_size / 10;

    let iter = sequence
        .windows(window_size)
        .step_by(step)
        .map(|window| {
            window
                .into_iter()
                .map(|c| match c {
                    b'G' | b'C' => 1,
                    _ => 0,
                })
                .sum::<u64>()
        })
        .map(move |sum| sum as f32 / window_size as f32);

    (0..).step_by(step).zip(iter)
}

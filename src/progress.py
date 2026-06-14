import time


def format_seconds_left(seconds_left):
    seconds_left = max(0, int(seconds_left))
    mins, secs = divmod(seconds_left, 60)
    if mins > 0:
        return f"~{mins} min {secs} s left"
    return f"~{secs} s left"


def update_method_progress(
    status_text,
    progress_bar,
    label,
    fraction,
    seconds_left,
    progress_start=0.0,
    progress_end=1.0,
):
    fraction = min(max(float(fraction), 0.0), 1.0)
    scaled_fraction = progress_start + (progress_end - progress_start) * fraction
    progress_bar.progress(scaled_fraction)
    status_text.info(
        f"{label}: **{fraction:.1%}** complete — {format_seconds_left(seconds_left)}"
    )


def make_estimated_progress_callback(status_text, progress_bar, progress_start, progress_end, label):
    def _callback(start_time, estimated_seconds):
        elapsed = time.time() - start_time
        fraction = min(elapsed / estimated_seconds, 0.95) if estimated_seconds else 0.95
        seconds_left = max(0, estimated_seconds - elapsed)
        update_method_progress(
            status_text,
            progress_bar,
            label,
            fraction,
            seconds_left,
            progress_start,
            progress_end,
        )

    return _callback


def make_done_total_progress_callback(status_text, progress_bar, progress_start, progress_end, label):
    def _callback(done, total, seconds_left):
        fraction = done / total if total else 1.0
        update_method_progress(
            status_text,
            progress_bar,
            label,
            fraction,
            seconds_left,
            progress_start,
            progress_end,
        )

    return _callback


def make_timed_done_total_progress_callback(status_text, progress_bar, progress_start, progress_end, label):
    def _callback(done, total, start_time):
        fraction = done / total if total else 1.0
        elapsed = time.time() - start_time
        estimated_total = elapsed / fraction if fraction > 0 else 0
        seconds_left = max(0, estimated_total - elapsed)
        update_method_progress(
            status_text,
            progress_bar,
            label,
            fraction,
            seconds_left,
            progress_start,
            progress_end,
        )

    return _callback

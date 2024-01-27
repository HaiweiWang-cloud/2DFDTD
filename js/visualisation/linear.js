function linear(max, min, value) {
    let val = Math.min(Math.max(value, min), max-0.0001);
    const d = max - min;
    val = d == 0 ? 0.5 : (val - min) / d;

    return [val*255, val*255, val*255, val*255];
}
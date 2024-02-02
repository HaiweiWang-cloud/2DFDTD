class GaussianSource extends Point {
    constructor(x, y, amplitude, duration, offset) {
        super(x, y);
        this.amplitude = amplitude;
        this.duration = duration;
        this.offset = offset;
    }

    getValue(t) {
        return this.amplitude * Math.exp(-1 * (t-this.offset)**2 / this.duration**2)
    }
}
class SineSource extends Point {
    constructor(x, y, amplitude, period, phase) {
        super(x, y);
        this.amplitude = amplitude;
        this.freq = 1/period;
        this.phase = phase;
    }

    getValue(t) {
        return this.amplitude * Math.sin(2*Math.PI*this.freq*t + this.phase);
    }
}
class SineSource extends Point {
    constructor(x, y, amplitude, period, phase) {
        super(x, y);
        this.amplitude = 0;
        this.maxAmplitude = amplitude;
        this.freq = 1/period;
        this.period = period;
        this.phase = phase;
        this.t = 0;
    }

    update(dt) {
        this.amplitude += 0.318 / this.period * this.maxAmplitude * Math.exp(-1 * ( (this.t-3*this.period) / (this.period) )**2);
        this.t += dt
    }

    getValue() {
        return this.amplitude * Math.sin(2*Math.PI*this.freq*this.t + this.phase);
    }
}
class Field1D {
    constructor(N, h, dt) {
        // Calculation of incident field\
        this.N = N;
        this.h = h;
        this.dt = dt;
        this.Cb = (this.dt/this.h)**2;
        this.hPML = 50;
        this.E = new Float32Array(N+this.hPML);
        this.H = new Float32Array(N+this.hPML);

        this.CaPML = []; // From interior to exterior
        this.CbPML = [];
        this.DaPML = [];
        this.DbPML = [];

        this.source = null;

        this.t = 0;

        this.calculatePML();
    }

    calculatePML() {
        for (let i=0; i < this.hPML; i++) {
            const sigma = ((i+1) / (this.hPML+1))**3;
            this.CaPML.push((1-0.5*sigma*this.dt) / (1+0.5*sigma*this.dt));
            this.CbPML.push((this.dt / this.h)**2 / (1+0.5*sigma*this.dt));
            this.DaPML.push((1-0.5*sigma*this.dt) / (1+0.5*sigma*this.dt));
            this.DbPML.push(1 / (1+0.5*sigma*this.dt));
        }
    }

    applySource() {
        this.E[0] = this.source.getValue(this.t);
    }

    update() {
        this.applySource();
        // H field update
        for (let i=0; i<this.N; i++) {
            this.H[i] += this.E[i+1] - this.E[i]; 
        }

        // PML region
        for (let t=0; t<this.hPML-1; t++) {
            this.H[this.N+t] = this.DaPML[t] * this.H[this.N+t] + this.DbPML[t] * (this.E[this.N+t+1] - this.E[this.N+t]);
        }

        // PEC boundary
        this.H[this.N+this.hPML-1] = this.DaPML[this.hPML-1] * this.H[this.N+this.hPML-1] + this.DbPML[this.hPML-1] * (0 - this.E[this.N+this.hPML-1]);

        // E field update
        for (let i=1; i<this.N; i++) {
            this.E[i] += this.Cb * (this.H[i] - this.H[i-1]);
        }

        // PML region

        for (let t=0; t<this.hPML; t++) {
            this.E[this.N+t] = this.CaPML[t] * this.E[this.N+t] + this.CbPML[t] * (this.H[this.N+t] - this.H[this.N+t-1]);
        }
        this.t++;
    }
}
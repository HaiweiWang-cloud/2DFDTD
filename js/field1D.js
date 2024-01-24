class Field1D {
    constructor(N, h, dt) {
        // Calculation of incident field\
        this.N = N;
        this.E = new Float32Array(N);
        this.H = new Float32Array(N);
        this.h = h;
        this.dt = dt;
        this.Cb = (this.dt/this.h)**2;
    }

    update() {
        for (let i=0; i<this.N; i++) {
            
        }
    }
}
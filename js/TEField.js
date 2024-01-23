class TMField {
    constructor(Nx, Ny, dt, h) {
        // Assuming nonpermeable media.
        this.Nx = Nx;
        this.Ny = Ny;
        this.dt = dt;
        this.h = h;
        this.num = Nx * Ny;

        this.Ez = new Float32Array(this.num); // Normalised
        this.Hx = new Float32Array(this.num);
        this.Hy = new Float32Array(this.num);
        this.mediaEz = new Uint8Array(this.num);

        this.Ca = [];
        this.Cb = [];

        this.n = 0;

        this.calculateUpdateCoefficients(1, 0);
    }

    calculateUpdateCoefficients(epsilonRZ, sigmaZ) {
        this.Ca.push((1-0.5*sigmaZ*this.dt/epsilonRZ)/(1+0.5*sigmaZ*this.dt/epsilonRZ));
        this.Cb.push((this.dt/this.h)**2/epsilonRZ/(1+0.5*sigmaZ*this.dt/epsilonRZ)); // Normalised
    }

    update() {
        const N = this.Ny;

        // Update H

        for (let i=0; i<this.Nx-1; i++) {
            for (let j=0; j<this.Ny-1; j++) {
                this.Hx[i*N+j] += this.Ez[i*N+j] - this.Ez[i*N+j+1];
                this.Hy[i*N+j] += this.Ez[(i+1)*N+j] - this.Ez[i*N+j];
            }
        }

        // Dirichlet boundary conditions
        // y-hi
        for (let i=0; i<this.Nx; i++) {
            this.Hx[i*N+N-1] += this.Ez[i*N+N-1] - 0;            
        }

        // x-hi
        for (let j=0; j<this.Ny; j++) {
            this.Hy[(this.Nx-1)*N+j] += 0 - this.Ez[(this.Nx-1)*N+j];
        }

        // Update Ez

        for (let i=1; i<this.Nx; i++) {
            for (let j=1; j<this.Ny; j++) {
                this.Ez[i*N+j] = this.Ca[this.mediaEz[i*N+j]] * this.Ez[i*N+j] + this.Cb[this.mediaEz[i*N+j]] * (this.Hy[i*N+j]-this.Hy[(i-1)*N+j]+this.Hx[i*N+j-1]-this.Hx[i*N+j]);
            }
        }

        // Dirichlet boundary conditions
        // y-lo
        for (let i=1; i<this.Nx; i++) {
            this.Ez[i*N] = this.Ca[this.mediaEz[i*N]] * this.Ez[i*N] + this.Cb[this.mediaEz[i*N]] * (this.Hy[i*N]-this.Hy[(i-1)*N]+ 0 -this.Hx[i*N]);
        }

        // x-lo
        for (let j=1; j<this.Ny; j++) {
            this.Ez[j] = this.Ca[this.mediaEz[j]] * this.Ez[j] + this.Cb[this.mediaEz[j]] * (this.Hy[j]- 0 +this.Hx[j-1]-this.Hx[j]);
        }

        //y-lo and x-lo
        this.Ez[0] = this.Ca[this.mediaEz[0]] * this.Ez[0] + this.Cb[this.mediaEz[0]] * (this.Hy[0] - 0 + 0 - this.Hx[0]);

        this.n++;
    }
}
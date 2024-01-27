class TMField {
    constructor(Nx, Ny, dt, h, periodicX, periodicY) {
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

        //this.incident = new Field1D(, this.h, this.dt)

        this.Ca = [];
        this.Cb = [];

        this.periodicX = periodicX ? 1 : 0;
        this.periodicY = periodicY ? 1 : 0;

        this.n = 0;

        this.calculateUpdateCoefficients(1, 0);
    }

    calculateUpdateCoefficients(epsilonRZ, sigmaZ) {
        this.Ca.push((1-0.5*sigmaZ*this.dt/epsilonRZ)/(1+0.5*sigmaZ*this.dt/epsilonRZ));
        this.Cb.push((this.dt/this.h)**2/epsilonRZ/(1+0.5*sigmaZ*this.dt/epsilonRZ)); // Normalised
    }

    setSFTF(tx0, tx1, ty0, ty1, phi) {
        this.tx0 = tx0;
        this.tx1 = tx1;
        this.ty0 = ty0;
        this.ty1 = ty1;
        this.phi = phi;
        this.cphi = Math.cos(phi);
        this.sphi = Math.sin(phi);
    }

    update() {
        const N = this.Ny; 
        const duration = 50;

        // Update H

        for (let i=0; i<this.Nx-1; i++) {
            for (let j=0; j<this.Ny-1; j++) {
                this.Hx[i*N+j] += this.Ez[i*N+j] - this.Ez[i*N+j+1];
                this.Hy[i*N+j] += this.Ez[(i+1)*N+j] - this.Ez[i*N+j];
            }
        }

        // Scattered field corrections
        /* for (let j=this.ty0-1; j<this.Ny-this.ty1; j++) {
            // x-lo, subtract incident Ez to right
            this.Hy[(this.tx0-1)*N+j] -= this.sampleIncidentE(this.tx0 * this.h, j * this.h);
            // x-hi, subtract incident Ez to left
            this.Hy[(this.Nx-this.tx1-1)*N+j] += this.sampleIncidentE((this.Nx-this.tx1-1) * this.h, j * this.h);
        } */
        
        for (let i=0; i<this.Nx-1; i++) {
            // y-lo
            this.Hx[i*N+this.ty0-1] += 1*Math.exp(-1 * (this.n-5*duration)**2 / (duration)**2);
            // y-hi
            //this.Hx[i*N+this.Ny-this.ty1-1] -= this.sampleIncidentE(i * this.h, (this.Ny-this.ty1-1) * this.h);
        }

        // Dirichlet boundary conditions
        // y-hi
        for (let i=0; i<this.Nx; i++) {
            const Eout = this.periodicY * this.Ez[i*N];
            this.Hx[i*N+N-1] += this.Ez[i*N+N-1] - Eout;            
        }

        // x-hi
        for (let j=0; j<this.Ny; j++) {
            const Eout = this.periodicX * this.Ez[j]
            this.Hy[(this.Nx-1)*N+j] += Eout - this.Ez[(this.Nx-1)*N+j];
        }

        // Update Ez

        for (let i=1; i<this.Nx; i++) {
            for (let j=1; j<this.Ny; j++) {
                this.Ez[i*N+j] = this.Ca[this.mediaEz[i*N+j]] * this.Ez[i*N+j] + this.Cb[this.mediaEz[i*N+j]] * (this.Hy[i*N+j]-this.Hy[(i-1)*N+j]+this.Hx[i*N+j-1]-this.Hx[i*N+j]);
            }
        }

        // Scattered field corrections
        /* for (let j=this.ty0-1; j<this.Ny-this.ty1; j++) {
            // x-lo, add incident to the left Hy field
            this.Ez[(this.tx0-1)*N+j] -= this.Cb[this.mediaEz[(this.tx0-1)*N+j]] * this.sampleIncidentH((this.tx0 - 1 -0.5) * this.h, j * this.h) * -this.cphi;
            // x-hi, add incident to the right Hy field
            this.Ez[(this.Nx - this.tx1 - 1)*N+j] += this.Cb[this.mediaEz[(this.Nx - this.tx1 - 1)*N+j]] * this.sampleIncidentH((this.Nx - this.tx1 - 1 + 0.5) * this.h, j * this.h) * -this.cphi;
        } */

        for (let i=0; i<this.Nx-1; i++) {
            // y-lo, add incident to the bottom Hx field
            this.Ez[i*N + this.ty0 - 1] += this.Cb[this.mediaEz[i*N + this.ty0 - 1]] * this.h / this.dt * Math.exp(-1 * (this.n-5*duration+this.h / this.dt*0.5-5)**2 / (duration)**2)// this.sampleIncidentH(i * this.h, (this.ty0 - 1 - 0.5) * this.h)* this.sphi;
            // y-hi, add incident to the top Hx field
            // this.Ez[i*N + this.Ny - this.ty1 - 1] -= this.Cb[this.mediaEz[i*N + this.Ny - this.ty1 - 1]] * this.sampleIncidentH(i * this.h, (this.Ny - this.ty1 - 1 + 0.5) * this.h) * this.sphi;
        }

        // Dirichlet boundary conditions
        // y-lo
        for (let i=1; i<this.Nx; i++) {
            const HxOut = this.periodicY * this.Hx[i*N+N-1];
            this.Ez[i*N] = this.Ca[this.mediaEz[i*N]] * this.Ez[i*N] + this.Cb[this.mediaEz[i*N]] * (this.Hy[i*N]-this.Hy[(i-1)*N] + HxOut -this.Hx[i*N]);
        }

        // x-lo
        for (let j=1; j<this.Ny; j++) {
            const HyOut = this.periodicX * this.Hy[(this.Nx-1) * N + j]
            this.Ez[j] = this.Ca[this.mediaEz[j]] * this.Ez[j] + this.Cb[this.mediaEz[j]] * (this.Hy[j] - HyOut +this.Hx[j-1]-this.Hx[j]);
        }

        //y-lo and x-lo
        this.Ez[0] = this.Ca[this.mediaEz[0]] * this.Ez[0] + this.Cb[this.mediaEz[0]] * (this.Hy[0] - this.periodicX * this.Hy[(this.Nx-1) * N] + this.periodicY * this.Hx[N-1] - this.Hx[0]);

        this.n++;
    }

    sampleIncidentH(x, y) {

    }

    sampleIncidentE(x, y) {

    }
}
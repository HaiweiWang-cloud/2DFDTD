class TMFieldPMLBerenger {
    constructor(Nx, Ny, dt, h) {
        // Assuming nonpermeable media.
        this.Nx = Nx;
        this.Ny = Ny;
        this.dt = dt;
        this.h = h;
        this.num = Nx * Ny;

        this.hPML = 20; // PML thickness

        this.Ezx = new Float32Array(this.num); // Normalised
        this.Ezy = new Float32Array(this.num);
        this.Hx = new Float32Array(this.num);
        this.Hy = new Float32Array(this.num);
        this.mediaEz = new Uint8Array(this.num);

        // Derived field quantities
        this.Ez = new Float32Array(this.num);
        this.intensity = new Float32Array(this.num);

        this.Ca = [];
        this.Cb = [];

        this.CaPML = []; // From interior to exterior
        this.CbPML = [];
        this.DaPML = [];
        this.DbPML = [];

        this.tfsfEnabled = false;

        this.n = 0;

        this.calculateUpdateCoefficients(1, 0);
        this.calculatePMLUpdateCoefficients();
    }

    findProjected(x, y) {
        return x * this.cphi + y * this.sphi;
    }

    sampleIncidentE(x, y) {
        const d = this.findProjected(x-this.x0 * this.h,y-this.y0 * this.h);
        return this.auxField.sampleE(d+2*this.h);
    }

    sampleIncidentH(x, y) {
        const d = this.findProjected(x-this.x0 * this.h, y-this.y0 * this.h);
        return this.auxField.sampleH(d+2*this.h);
    }

    updateDerived() {
        this.Ezx.forEach((xComp, index) => {
            this.Ez[index] = xComp + this.Ezy[index];
        });
        this.Ez.forEach((E, index) => {
            this.intensity[index] = E * E * this.Ca[this.mediaEz[index]];
        });
    }

    calculateUpdateCoefficients(epsilonRZ, sigmaZ) {
        this.Ca.push((1-0.5*sigmaZ*this.dt/epsilonRZ)/(1+0.5*sigmaZ*this.dt/epsilonRZ));
        this.Cb.push((this.dt/this.h)**2/epsilonRZ/(1+0.5*sigmaZ*this.dt/epsilonRZ)); // Normalised
    }

    calculatePMLUpdateCoefficients() {
        for (let i=0; i < this.hPML; i++) {
            const sigma = ((i+1) / (this.hPML+1))**3;
            this.CaPML.push((1-0.5*sigma*this.dt) / (1+0.5*sigma*this.dt));
            this.CbPML.push((this.dt / this.h)**2 / (1+0.5*sigma*this.dt));
            this.DaPML.push((1-0.5*sigma*this.dt) / (1+0.5*sigma*this.dt));
            this.DbPML.push(1 / (1+0.5*sigma*this.dt));
        }
    }

    setTFSF(phi, x0, x1, y0, y1) {
        this.auxField = new Field1D(2+this.Nx+this.Ny+50, this.h, this.dt);
        this.phi = phi;
        this.cphi = Math.cos(phi);
        this.sphi = Math.sin(phi);
        this.auxField.source = new SineSource(N-1, 0, 1, 600, 0);
        this.tfsfEnabled = true;
        // boundaries of the total-field rectangle
        this.x0 = x0;
        this.y0 = y0;
        this.x1 = x1;
        this.y1 = y1;
    }

    computeTFSFE() {
        // x-direction
        for (let j=this.y0; j<this.y1+1; j++) {
            // x-lo, add incident to the left Hy field
            this.Ezx[this.x0*N+j] -= this.Cb[this.mediaEz[this.x0*N+j]] * this.sampleIncidentH((this.x0 - 0.5) * this.h, j * this.h) * -this.cphi;
            
            // x-hi, add incident to the right Hy field
            this.Ezx[this.x1*N+j] += this.Cb[this.mediaEz[this.x1*N+j]] * this.sampleIncidentH((this.x1 + 0.5) * this.h, j * this.h) * -this.cphi;
            
        }
    }

    computeTFSFH() {
        // x-direction
        for (let j=this.y0; j<this.y1+1; j++) {
            // x-lo, subtract incident Ez to right
            this.Hy[(this.x0-1)*N+j] -= this.sampleIncidentE(this.x0 * this.h, j * this.h);
            // x-hi, subtract incident Ez to left
            this.Hy[this.x1*N+j] += this.sampleIncidentE(this.x1 * this.h, j * this.h);
        }   
    }

    update() {
        const N = this.Ny; 
        const duration = 50;

        // Update H physical region

        for (let i=this.hPML-1; i<this.Nx-this.hPML; i++) {
            for (let j=this.hPML-1; j<this.Ny-this.hPML; j++) {
                this.Hx[i*N+j] += this.Ezx[i*N+j] - this.Ezx[i*N+j+1] + this.Ezy[i*N+j] - this.Ezy[i*N+j+1];
                this.Hy[i*N+j] += this.Ezx[(i+1)*N+j] - this.Ezx[i*N+j] + this.Ezy[(i+1)*N+j] - this.Ezy[i*N+j];
            }
        }

        // Update H PML
        // y-direction
        for (let i=0; i<this.Nx-1; i++) {
            // y-low
            for (let t=0; t<this.hPML; t++) {
                const Da = this.DaPML[t];
                const Db = this.DbPML[t];
                this.Hx[i*N+this.hPML-t-1] = Da * this.Hx[i*N+this.hPML-t-1] + Db * (this.Ezx[i*N+this.hPML-t-1] - this.Ezx[i*N+this.hPML-t] + this.Ezy[i*N+this.hPML-t-1] - this.Ezy[i*N+this.hPML-t]);
            }
            
            // y-hi
            for (let t=0; t<this.hPML-1; t++) {
                const Da = this.DaPML[t];
                const Db = this.DbPML[t];
                this.Hx[i*N+(this.Ny-this.hPML+t)] = Da * this.Hx[i*N+(this.Ny-this.hPML+t)] + Db * (this.Ezx[i*N+(this.Ny-this.hPML+t)] - this.Ezx[i*N+(this.Ny-this.hPML+t)+1] + this.Ezy[i*N+(this.Ny-this.hPML+t)] - this.Ezy[i*N+(this.Ny-this.hPML+t)+1])
            }

            // y-hi PEC boundary
            this.Hx[i*N+N-1] = this.DaPML[this.hPML-1] * this.Hx[i*N+N-1] + this.DbPML[this.hPML-1] * (this.Ezx[i*N+N-1] + this.Ezy[i*N+N-1] - 0);
        }
            

        // x-direction
        for (let j=0; j<this.Ny-1; j++) {
            // x-low
            for (let t=0; t<this.hPML; t++) {
                const Da = this.DaPML[t];
                const Db = this.DbPML[t];
                this.Hy[(this.hPML-t-1)*N+j] = Da * this.Hy[(this.hPML-t-1)*N+j] + Db * (this.Ezx[(this.hPML-t)*N+j] - this.Ezx[(this.hPML-t-1)*N+j] + this.Ezy[(this.hPML-t)*N+j] - this.Ezy[(this.hPML-t-1)*N+j]);
            }

            // x-hi
            for (let t=0; t<this.hPML-1; t++) {
                const Da = this.DaPML[t];
                const Db = this.DbPML[t];
                this.Hy[(this.Nx-this.hPML+t)*N+j] = Da * this.Hy[((this.Nx-this.hPML+t))*N+j] + Db * (this.Ezx[(this.Nx-this.hPML+t+1)*N+j] - this.Ezx[(this.Nx-this.hPML+t)*N+j] + this.Ezy[(this.Nx-this.hPML+t+1)*N+j] - this.Ezy[(this.Nx-this.hPML+t)*N+j]);
            }

            // x-hi PEC boundary
            this.Hy[(this.Nx-1)*N+j] = this.DaPML[this.hPML-1] * this.Hy[(this.Nx-1)*N+j] + this.DbPML[this.hPML-1] * ( 0 - this.Ezx[(this.Nx-1)*N+j] - this.Ezy[(this.Nx-1)*N+j] );
        }

        // Update E Physical
        for (let i=this.hPML-1; i<this.Nx-this.hPML; i++) {
            for (let j=this.hPML-1; j<this.Ny-this.hPML; j++) {
                this.Ezx[i*N+j] = this.Ca[this.mediaEz[i*N+j]] * this.Ezx[i*N+j] + this.Cb[this.mediaEz[i*N+j]] * (this.Hy[i*N+j]-this.Hy[(i-1)*N+j]);
                this.Ezy[i*N+j] = this.Ca[this.mediaEz[i*N+j]] * this.Ezy[i*N+j] + this.Cb[this.mediaEz[i*N+j]] * (this.Hx[i*N+j-1]-this.Hx[i*N+j]);
            }
        }
        
        // Update E PML
        // y direction
        for (let i=1; i<this.Nx; i++) {
            // y-low
            for (let t=1; t<this.hPML; t++) {
                const Ca = this.CaPML[t];
                const Cb = this.CbPML[t];
                this.Ezy[i*N+this.hPML-t] = Ca * this.Ezy[i*N+this.hPML-t] + Cb * (this.Hx[i*N+this.hPML-t-1]-this.Hx[i*N+this.hPML-t]);
            }

            // y-low PMC boundary
            this.Ezy[i*N] = this.CaPML[this.hPML-1] * this.Ezy[i*N] + this.CbPML[this.hPML-1] * (0 -this.Hx[i*N]);

            // y-hi
            for (let t=0; t<this.hPML; t++) {
                const Ca = this.CaPML[t];
                const Cb = this.CbPML[t];
                this.Ezy[i*N+this.Ny-this.hPML+t] = Ca * this.Ezy[i*N+this.Ny-this.hPML+t] + Cb * (this.Hx[i*N+this.Ny-this.hPML+t-1]-this.Hx[i*N+this.Ny-this.hPML+t]);
            }
        }

        // x direction
        for (let j=1; j<this.Ny; j++) {
            // x-low
            for (let t=1; t<this.hPML; t++) {
                const Ca = this.CaPML[t];
                const Cb = this.CbPML[t];
                this.Ezx[(this.hPML-t)*N+j] = Ca * this.Ezx[(this.hPML-t)*N+j] + Cb * (this.Hy[(this.hPML-t)*N+j]-this.Hy[(this.hPML-t-1)*N+j]);
            }

            // x-low PMC boundary
            this.Ezx[j] = this.CaPML[this.hPML-1] * this.Ezx[j] + this.CbPML[this.hPML-1] * (this.Hy[j] - 0);

            // x-hi
            for (let t=0; t<this.hPML; t++) {
                const Ca = this.CaPML[t];
                const Cb = this.CbPML[t];
                this.Ezx[(this.Nx-this.hPML+t)*N+j] = Ca * this.Ezx[(this.Nx-this.hPML+t)*N+j] + Cb * (this.Hy[(this.Nx-this.hPML+t)*N+j]-this.Hy[(this.Nx-this.hPML+t-1)*N+j]);
            }
        }

        this.n++;
    }
}

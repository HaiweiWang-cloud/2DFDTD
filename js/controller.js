const TE = 0;
const TM = 1;

class Controller {
    constructor(canvas) {
        // Field
        this.sources = [];
        
        // Display
        this.canvas = canvas;
        this.ctx = canvas.getContext("2d");
        this.cmScalar = new Twilight();
        this.cmIntensity = new Magma();

        this.maxAmplitude = 1;

        // Editor
        this.mouse = null;
        this.hovered = null;
        this.selected = null;
        this.dragging = false;

        this.#addEventListeners();
    }

    initialiseField(polarisation, Nx, Ny, dt, h) {
        if (polarisation == TE) {
            this.field = new TEField(Nx, Ny, dt, h);
        } else if (polarisation == TM) {
            this.field = new TMFieldPMLBerenger(Nx, Ny, dt, h);
        }

        this.pxsX = Math.floor(this.canvas.width / Nx);
        this.pxsY = Math.floor(this.canvas.height / Ny);
    }

    applySources() {
        const N = this.field.Ny;
        for (const source of this.sources) {
            source.update(1);
            this.field.Ezx[Math.floor(source.x*N+source.y)] += 0.5*source.getValue(this.field.n);
            this.field.Ezy[Math.floor(source.x*N+source.y)] += 0.5*source.getValue(this.field.n);
        }
    }

    draw({intensity=false, vector=false} = {}) {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        this.field.updateDerived();
        if (intensity) {
            this.drawIntensity();
        } else {
            this.drawScalar();
        }

        if (vector) {
            this.drawVector();
        }

        this.drawSelectedPoints();
    }

    drawScalar() {
        colormesh(this.ctx, this.canvas.width, this.canvas.height, this.field.Ez, this.field.Nx, this.field.Ny, this.cmScalar.getValue.bind(this.cmScalar), {max: this.maxAmplitude, min: -this.maxAmplitude});
    }

    drawIntensity() {
        colormesh(this.ctx, this.canvas.width, this.canvas.height, this.field.intensity, this.field.Nx, this.field.Ny, this.cmIntensity.getValue.bind(this.cmIntensity), {min: 0, max: this.maxAmplitude * this.maxAmplitude});
    }

    drawVector() {
        vector2D(this.ctx, this.canvas.width, this.canvas.height, this.field.Hx, this.field.Hy, this.field.Nx, this.field.Ny, 16, 16, 3, {lineWidth: 1});
    }

    drawSelectedPoints() {
        if (this.selected) {
            const point = new Point(this.selected.x * this.pxsX, this.selected.y * this.pxsY);
            point.draw(this.ctx, {size: 5, color: "white", outline: true});
        }
        if (this.hovered) {
            const point = new Point(this.hovered.x * this.pxsX, this.hovered.y * this.pxsY);
            point.draw(this.ctx, {size: 5, color: "purple", outline: true});
        }
    }

    #addEventListeners() {
        this.canvas.addEventListener("mousemove", this.#handleMouseMove.bind(this));
        this.canvas.addEventListener("mousedown", this.#handleMouseDown.bind(this));
        this.canvas.addEventListener("mouseup", () => this.dragging = false);
        this.canvas.addEventListener("touchend", () => this.dragging = false);
        this.canvas.addEventListener("contextmenu", (evt) => evt.preventDefault());
    }

    #handleDrag() {
        if (this.dragging && this.selected) {
            this.selected.x = this.mouse.x;
            this.selected.y = this.mouse.y;
        }
    }

    #handleMouseMove(evt) {
        this.mouse = new Point(evt.offsetX / this.pxsX, evt.offsetY / this.pxsY);
        this.hovered = getNearestPoint(this.mouse, this.sources, 10);
    
        this.#handleDrag();
    }

    #handleMouseDown(evt) {
        this.dragging = true;
        if (evt.button == 0) {
            this.#handleLeftClick();
        }

        if (evt.button == 2) {
            if (this.selected) {
                this.selected = null;
            }
        }
    }

    #handleLeftClick() {
        if (this.hovered) {
            this.selected = this.hovered;
        } else if (this.selected) {
            this.selected = null;
        } else {
            this.selected = new SineSource(this.mouse.x, this.mouse.y, 4, 300, 0);
            this.sources.push(this.selected);
        }
    }
}
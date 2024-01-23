const TE = 0;
const TM = 1;

class Controller {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext("2d");
        this.cmScalar = new Twilight();

        this.sources = [];

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
            this.field = new TMField(Nx, Ny, dt, h);
        }

        this.pxsX = Math.floor(this.canvas.width / Nx);
        this.pxsY = Math.floor(this.canvas.height / Ny);
    }

    applySources() {
        const N = this.field.Ny;
        for (const source of this.sources) {
            this.field.Ez[Math.floor(source.x*N+source.y)] = source.getValue(this.field.n);
        }
    }

    draw() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        this.drawScalar();
        this.drawSelectedPoints();
    }

    drawScalar() {
        colormesh(this.ctx, this.canvas.width, this.canvas.height, this.field.Ez, this.field.Nx, this.field.Ny, this.cmScalar.getValue.bind(this.cmScalar), {max: 1, min: -1});
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
            this.selected = new GaussianSource(this.mouse.x, this.mouse.y, 2, 20, this.field.n + 6*20);
            this.sources.push(this.selected);
        }
    }
}
class Controller1D {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext("2d");
        this.field = null;
    }

    initialiseField(N, h, dt) {
        this.field = new Field1D(N, h, dt);
    }

    setSource(N) {
        this.field.source = new SineSource(N-1, 0, 1, 100, 0);
    }

    draw() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        lineplot(this.ctx, this.field.E, this.canvas.width, this.canvas.height, {vmin: -1, vmax: 1});
    }
}
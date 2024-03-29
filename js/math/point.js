class Point {
    constructor(x, y) {
        this.x = x;
        this.y = y;
    }

    equals(point) {
        return (this.x === point.x && this.y === point.y)
    }

    draw(ctx, { size = 18, color = "black", outline = false, fill = false} = {}) {
        const rad = size / 2;
        if (outline) {
            ctx.beginPath();
            ctx.strokeStyle = color;
            ctx.arc(this.x, this.y, rad, 0, 2 * Math.PI);
            ctx.stroke();
        }
        if (fill) {
            ctx.beginPath();
            ctx.arc(this.x, this.y, rad, 0, 2 * Math.PI);
            ctx.fillStyle = color;
            ctx.fill();
        }
    }
}
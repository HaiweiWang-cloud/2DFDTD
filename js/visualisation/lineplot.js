function lineplot(ctx, yVals, Nx, Ny, {vmin=0, vmax=1, color="white", linewidth=3}={}) {
    const numPoints = yVals.length;
    const xScale = Nx / numPoints;
    let d = vmax - vmin;
    d = d==0 ? 0 : Ny / d;
   
    ctx.beginPath();
    ctx.strokeStyle = color;
    ctx.lineWidth = linewidth;
    ctx.moveTo(0, Math.min(vmax, Math.max(Ny, Math.min(0, (yVals[0]-vmin) * d))));
    yVals.forEach((y, i) => {
        ctx.lineTo(i*xScale, Math.min(Ny, Math.max(0, (y-vmin) * d)));
    });
    ctx.stroke();
}
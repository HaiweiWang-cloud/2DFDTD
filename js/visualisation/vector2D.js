function vector2D (ctx, width, height, vx, vy, Nx, Ny, skipX, skipY, scale, {color="white", lineWidth=2}={}) {
    ctx.strokeStyle = color;
    ctx.lineWidth = lineWidth;
    const lx = Math.floor(width / Nx);
    const ly = Math.floor(height / Ny);

    for (let x = skipX-1; x<Nx-skipX; x+=skipX) {
        for (let y= skipY-1; y<Ny-skipY; y+=skipY) {
            const x0 = x*lx;
            const y0 = y*ly;
            const vxVal = vx[x*Ny+y];
            const vyVal = vy[x*Ny+y];
            const nx = skipX * lx * vxVal * scale;
            const ny = skipY * ly * vyVal * scale;
            ctx.beginPath();
            ctx.moveTo(x0, y0);
            ctx.lineTo(x0+nx, y0+ny);
            ctx.arc(x0+nx, y0+ny, lineWidth*0.6, 0, 2*Math.PI);
            ctx.stroke();
        }
    }
}
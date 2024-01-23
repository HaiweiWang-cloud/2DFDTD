function colormesh(ctx, width, height, field, Nx, Ny, colormap, {alpha=1, min=null, max=null}={}) {
    const pixelSize = Math.floor(Math.min(width/Nx, height/Ny));
    let minF = 0;
    let maxF = 0;
    if (!min) {
        field.forEach((num) => {
            minF = Math.min(num, minF);
            
        });
        
    } else {
        minF = min;
    }

    if (!max) {
        field.forEach((num) => {
            maxF = Math.max(num, maxF);
        });
    } else {
        maxF = max;
    }

    const id = ctx.getImageData(0, 0, width, height);

    for (let i=0; i<Nx; i++) {
        for (let j=0; j<Ny; j++) {
            const color = colormap(maxF, minF, field[i*Ny + j]);
            
            const r = color[0];
            const g = color[1];
            const b = color[2];
            const a = color[3];

            const x = Math.floor(i * pixelSize);
            const y = Math.floor(j * pixelSize);
            for (let yi = y; yi < y + pixelSize; yi++) {
                let p = 4 * (yi * width + x);
                for (let xi = x; xi < x + pixelSize; xi++){
                    id.data[p++] = r;
                    id.data[p++] = g;
                    id.data[p++] = b;
                    id.data[p++] = a * alpha;
                }
            }
        }
    }

    ctx.putImageData(id, 0, 0);
    
}
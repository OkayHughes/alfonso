function cube = cubeFromBasis(basis)
    cube.n = basis.n;
    cube.d = basis.d/2;
    cube.L = basis.L;
    cube.U = basis.U;
    cube.w = basis.w;
    cube.pts = basis.pts;
    cube.P0 = basis.P0;
    [cube.P,~] = qr(basis.P0,0);
end
function [solP, matP, rhsP] = computeSolProj2D_DG(mesh, dofm)

matP = buildMatrixGlo2D_DG(mesh, dofm);
rhsP = buildVectorGloRhs2D_DG(mesh, dofm, @mySol);
solP = matP\rhsP;

end
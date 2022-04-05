function [solP, matP, rhsP] = computeSolProj2D_CG(mesh, dofm)

matP = buildMatrixGlo2D_CG(mesh, dofm);
rhsP = buildVectorGloRhs2D_CG(mesh, dofm, @mySol);
solP = matP\rhsP;

end
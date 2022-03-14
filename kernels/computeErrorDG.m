function [errorL2, errorH1] = computeErrorDG(mesh, dofm, solNum, solRef)

errNum = solNum(1:dofm.numDofTRI)-solRef;

[matM, matK] = buildMatrixGloDG(mesh, dofm);

normL2 = real(sqrt(solRef'*matM*solRef));
normH1 = real(sqrt(solRef'*matK*solRef + solRef'*matM*solRef));

errorL2 = real(sqrt(errNum'*matM*errNum)) / normL2;
errorH1 = real(sqrt(errNum'*matK*errNum + errNum'*matM*errNum)) / normH1;

end
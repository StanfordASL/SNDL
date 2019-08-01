function out = dF_dtheta_p_h_tom(alpha,lambda,x,tempD)
% dF_dtheta_p_h_tom - Autogenerated file.
%
% out = dF_dtheta_p_h_tom(alpha,lambda,x,tempD)
tempC4   = x((1:4)')'*tempD{1};
tempC5   = sin(tempC4);
tempC6   = cos(tempC4);
clear tempC4
tempC9   = x'*tempD{2};
tempC11  = kron(sin(tempC9),[0 1]);
tempC13  = kron(cos(tempC9),[1 0]);
clear tempC9
tempC15  = ekron(tempC13+tempC11,6,2,'I','I');
clear tempC11 tempC13
tempC26  = (((0.1443375672974065*tempC15(4,:))*alpha)*0.1666666666666667)*(kron(-(tempC5.*tempD{3}),[1 0])+kron(tempC6.*tempD{3},[0 1]));
tempC27  = ekron(tempC26,10,1,'I','I');
tempC39  = (((0.1443375672974065*tempC15(3,:))*alpha)*0.1666666666666667)*(kron(-(tempC5.*tempD{4}),[1 0])+kron(tempC6.*tempD{4},[0 1]));
clear tempC15
tempC40  = ekron(tempC39,10,1,'I','I');
tempC63  = tempD{8}*x;
tempC66  = -scalerows(sin(tempC63),tempD{8});
tempC69  = scalerows(cos(tempC63),tempD{8});
clear tempC63
tempC72  = 0.1443375672974065*(kron(tempC66,[1;0])+kron(tempC69,[0;1]));
clear tempC66 tempC69
tempC74  = alpha';
tempC86  = [vec(tempC74*kron(tempC72,tempD{12}));vec(tempC74*kron(tempC72,tempD{11}));vec(tempC74*kron(tempC72,tempD{10}));vec(tempC74*kron(tempC72,tempD{9}))];
clear tempC72 tempC74
tempC90  = kron(tempC6,[1 0])+kron(tempC5,[0 1]);
clear tempC5 tempC6
tempC91  = 0.1666666666666667*tempC90;
clear tempC90
tempC92  = ekron(tempC91,10,1,'I','I');
tempC93  = setrows(21,(1:10)',tempC92);
tempC94  = tempC93((7:9)',:);
tempC95  = tempC93([4;5],:);
tempC96  = tempC92((7:10)',:);
tempC97  = tempC92((4:6)',:);
tempC98  = tempC92([2;3],:);
clear tempC92
tempC99  = tempC91*tempD{5};
tempC100 = tempC91*tempD{6};
clear tempC91
out      = ((ekron(sparse(tempD{14},tempD{13},tempC86,4,6),4,1,'I',submatrix(setrows(27,1,tempC100,[2;3],tempC98,(4:6)',tempC97,(7:10)',tempC96,22,tempC99,[23;24],tempC95,(25:27)',tempC94),tempD{15},':'))-submatrix(setrows(16,1,tempC39*tempD{6}+tempC26*tempD{6},[2;3],tempC40([2;3],:)+tempC27([2;3],:),(4:6)',tempC40((4:6)',:)+tempC27((4:6)',:),(7:10)',tempC40((7:10)',:)+tempC27((7:10)',:),11,tempC39*tempD{5}+tempC26*tempD{5},[12;13],tempC40([4;5],:)+tempC27([4;5],:),(14:16)',tempC40((7:9)',:)+tempC27((7:9)',:)),tempD{7},':'))+reshape((reshape(tempC86,6,4)'*reshape(reshape(submatrix(setrows(25,1,tempC100,[2;3],tempC98,(4:6)',tempC97,(7:10)',tempC96,11,tempC99,[12;13],tempC95,(14:16)',tempC94,(17:20)',tempC93((11:14)',:),(21:25)',tempC93((16:20)',:)),tempD{16},':'),4,4320)',6,2880))',720,16)')+(2*(lambda+0.1))*submatrix(setrows(16,1,tempC100,[2;3],tempC98,(4:6)',tempC97,(7:10)',tempC96,11,tempC99,[12;13],tempC95,(14:16)',tempC94),tempD{7},':');
function Q = Qobj3(D, b)

over = D.initX1;

over(1)  = b(1);
over(2)  = b(2);
over(3)  = b(3);
over(5)  = b(4);
over(7)  = b(5);
over(8)  = b(6);
over(9)  = b(7);
over(11) = b(8);

[Smom,~,~] = modelF(over,D.fid);
diff = (D.Dmom-Smom);
diff = diff(D.moms);
Q = sqrt(diff'*D.W*diff);

if (~isnan(Q))
   % Print params being used
   fprintf(D.fid,'\nParams:\n');
   for i=1:numel(over)
      fprintf(D.fid,'%7.3f\n', over(i));
   end
   % print resulting Q
   fprintf(D.fid,'\nQ = %12.3f\n\n',Q);
end

end

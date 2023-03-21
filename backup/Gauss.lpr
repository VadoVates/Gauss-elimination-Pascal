program Gauss;

type
  tMacierz = array of array of double;

function EliminacjaGaussa (macierz : tMacierz; ileNiewiadomych : word; eps : double) : boolean;
var
  i, j, k : word;
  x1, x2 : double;
begin
  for i:=0 to ileNiewiadomych-2 do
  begin
    for j:=i+1 to ileNiewiadomych-1 do
    begin
      if (macierz[i][i] < eps) then
      begin
        EliminacjaGaussa := false;
        Exit;
      end;

    end;
  end;
end;

procedure Implementacja;
var
  macierzAB : tMacierz;
  ileNiewiadomych, i, j : word;
  eps : double;
begin
  write ('Podaj ile chcesz niewiadomych: ');
  read (ileNiewiadomych);

  write ('Podaj jakiej oczekujesz dokladnosci: ');
  read (eps);

  SetLength (macierzAB, ileNiewiadomych, ileNiewiadomych+1);
  for i:=0 to ileNiewiadomych-1 do
  begin
    for j:=0 to ileNiewiadomych do
    read (macierzAB[i][j]);
  end;
  for i:=0 to ileNiewiadomych-1 do
  begin
    for j:=0 to ileNiewiadomych do
    write (macierzAB[i][j]:3:0);
    writeln;
  end;

  if (EliminacjaGaussa (macierzAB, ileNiewiadomych, eps)) then
  begin

  end;
  else
  begin
    writeln ('Macierz osobliwa, det = 0');
  end;
end;

begin
  Implementacja;
  readln;
  readln;
end.

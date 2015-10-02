int ispow2(n)
int n;

{
  if(n < 2)
    return(0);
  for(;!(n%2);n /= 2);
  return(((n == 1)?(1):(0)));
}

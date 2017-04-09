program tsmr_main
use tsmr
implicit none
integer :: ipar(4),n,u,na,np,i, k_planets
REAL :: rpar(2),rtol,atol
Character(len=1000) :: fpar(7),cfg
logical :: ln

call get_command_argument(1,cfg)
call tsmr_cfg(cfg,ipar,rpar,ln,fpar)
n=ipar(1)
u=ipar(2)
na=ipar(3)
np=ipar(4)
rtol=rpar(1)
atol=rpar(2)
k_planets = int((-11+int(sqrt(real(169+8*n))))/2)

call tsmr_integ(n,u,na,np,ln,rtol,atol,fpar, k_planets)
end program tsmr_main
PROGRAM RaicesComplejas

IMPLICIT NONE

COMPLEX(8) :: x0
REAL(8),PARAMETER :: tolerancia= 0.0001
INTEGER,PARAMETER :: maxiter=1000
  
  x0=(1,1) !cargar valor de x0
                           
OPEN(UNIT=1,FILE='resolucion.dat',STATUS='REPLACE')
   CALL newtonraphson(x0)
   WRITE(*,*)'La resolucion se encuentra en archivo resolucion.dat'
CLOSE(UNIT=1,STATUS='KEEP')

CONTAINS

FUNCTION f(x) !cargar la Funcion

COMPLEX(8) x,f
   f=x**5-3.5*x**4+2.75*x**3+2.125*x**2-3.875*x+1.25

END FUNCTION f

FUNCTION fprima(x) !Cargar la derivada de la funcion
COMPLEX(8) x,fprima
    fprima=  5*x**4-14*x**3+8.25*x**2+4.25*x-3.875
END FUNCTION fprima

SUBROUTINE newtonraphson (x0)
  
  REAL(8) ::error, errorenx
  INTEGER :: iter
  COMPLEX(8) x0,x,xant

  WRITE(1,'(A/)') ' RESOLUCION MEDIANTE EL METODO DE NEWTONRAPSHON '
  error=(tolerancia*2.0)
  errorenx=(tolerancia*2.0)
  iter=1
  WRITE(*,*)'El valor de "x0" es ',x0
  
  
  x=x0
  DO WHILE (error>=tolerancia)  
!   DO WHILE (errorenx> tolerancia)
    WRITE(1,*)'---------------------------------------'
    WRITE(1,*)
    xant=x
    x=x-(f(x)/fprima(x))
    WRITE(1,'(A,I4)')'Iteracion: ',iter
    WRITE(1,'(A,F20.7)')'Valor de x: ',x
    error=abs(f(x))
    WRITE(1,'(A,F20.7)')'Error o distancia a y: ',error
    iter=iter+1
    errorenx= abs(x-xant)
    WRITE(1,'(A,F20.7)')'Error en x : ',errorenx

  END DO
  WRITE(1,*)'---------------------------------------'
  WRITE(1,*)
  WRITE(1,'(A,F20.7/)')'La raiz encontrada es: ',x
  WRITE(*,*)'La raiz encontrada es: ',x

END SUBROUTINE newtonraphson

END PROGRAM

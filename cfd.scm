;; test functions to evolve
(include "/home/dsp/devel/mhd-lisp/sagdeev.scm")
(define genstepfunction 
        (lambda (x1 x2 ymax ybase)
                (cond ((>= x1 x2) 'ERROR)
                      (else (lambda (x)
                                    (cond ((and (>= x x1) (<= x x2)) ymax)
                                          (else ybase)))))))

(define teststep (genstepfunction 0.5 1 2 1))
(define testxspace (genpoints 0 10 85))
(define testpulse (map teststep testxspace))
(define testsinpulse (map sin testxspace))

(define (nderiv lst tmp dx)
  (if (eq? (cddr lst) '())
      tmp
      (nderiv (cdr lst) (cons (list (caar lst) (/ (- (cadadr lst) (cadar lst)) (square dx))) tmp) dx)))
    

(define list-set! (lambda (lst k val) (set-car! (list-tail lst k) val)))'

(define (lst2gpstring lst)
  (let ((res ""))
    (cond 
     ((pair? (car lst))
        (begin
        (do-times n (length lst)
          (set! res
                (string-append res
                               (string-append
                                 (number->string (car (list-ref lst n))) "\t"
                                 (number->string (cadr (list-ref lst n))) "\n"))))
        (set! res (string-append res " e \n"))))
     (else
        (begin
        (do-times n (length lst)
         (set! res 
               (string-append res
                              (string-append 
                                (number->string n) "\t"
                                (number->string (list-ref lst n)) "\n"))))
        (set! res (string-append res " e \n")))))
  res))

               
          

;; step1: linear convection with constant velocity c. this is linear
;; and we solve using forward scheme for time and backwards for space
(define linearConv (lambda (pulse nt)
  (let* ((c 1)
         (sigma 0.5) ;; courant number for CFL
         (nx (length pulse))
         (dx (/ 2 (- nx  1)))
         (dt (* sigma dx))
	 (outpulse (list-copy pulse)))
    (do-times t nt
      (let ((copy (list-copy outpulse)))
        (do-times x nx
          (if (> x 0);ignore first time
            (list-set! outpulse x (- (list-ref copy x) (* c (/ dt dx) (- (list-ref copy x) (list-ref copy (- x 1))))))))))
    outpulse)))

(define nonlinearConv (lambda (pulse nt)
  (let* ((sigma 0.5)
         (nx (length pulse))
         (dx (/ 2 (- nx  1)))
         (dt (* sigma dx))
	 (outpulse (list-copy pulse)))
    (do-times t nt
      (let ((copy (list-copy outpulse)))
        (do-times x nx
          (if (> x 0);ignore first time
            (list-set! outpulse x (- (list-ref copy x) (* (list-ref copy x) (/ dt dx) (- (list-ref copy x) (list-ref copy (- x 1))))))))))
    outpulse)))

;; helper to only get the vals (y) of funs
(define (extractvals lst) (map (lambda (x) (cadr x)) lst))

;; how i got the velocities
;; (define vroots (map (lambda (x) (list (car x) (- 1.1 (sqrt (- (square 1.1) (* 2 (cadr x)))))))  sortedphiroots))
;; how i got the densities
;; (define nroots (map (lambda (x) (list (car x) (/ 1 (sqrt (- 1 (/ (* 2 (cadr x)) (square 1.1))))))) sortedphiroots))
;; how i sorted
;; (sort! phiroots (lambda (a b) (< (car a) (car b)))

(define (square x) (* x x))
;; takes argumens of n density u velocity f phi and the number of timesteps
(define plasmasystem (lambda (n u f nt)
  (if (not (and (eq? (length n) (length u)) (eq? (length n) (length f))))
      'ERRLISTFORM
      (let* ((sigma 0.5)
            (nx (length n))
            (dx (/ 2 (- nx 1)))
            (dt (* sigma dx))
            (n1 (list-copy n))
            (u1 (list-copy u))
            (f1 (list-copy f)))
        (do-times t nt
          (let ((ntemp (list-copy n1))
                (utemp (list-copy u1))
                (ftemp (list-copy f1)))
            (format #t "doing time:~A~%" t)
            (do-times x nx
              (format #t " point:~A " x)
              (if (< x (- nx 1)) ;; cause we're doing forward. ignore the last elem
                  (list-set! n1 x (- (list-ref ntemp x) (* (/ dt dx) (- (* (list-ref ntemp (+ x 1)) (list-ref utemp (+ x 1))) (* (list-ref ntemp x) (list-ref utemp x)))))))
              (if (< x (- nx 1))
                  (list-set! u1 x (- (list-ref utemp x) (* (/ dt dx) (+ (list-ref ftemp (+ x 1)) (list-ref ftemp x) (square (list-ref utemp (+ x 1))) (- (list-ref utemp x)))))))
              (if (and (> x 0) (< x (- nx 1))) ;; central in space
                  (list-set! f1 x (/ (- (* (square dx) (- 1 (list-ref ntemp x))) (+ (list-ref ftemp (+ x 1)) (list-ref ftemp (- x 1)))) (- (+ 2 (square dx)))))))))
        (list n1 u1 f1)))))

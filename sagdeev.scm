(include "mathh-constants")
(use format)
(use loops)
(use srfi-1)
(use srfi-13)
(use srfi-69)

(define sagdeev 
	(lambda (phi M)
		(+ 1 (- (expt E phi)) (* (* M M) (- 1 (sqrt (- 1 (/ (* 2 phi) (* M M)))))))))

(define DsagdeevDphi
	(lambda (phi M)
		(- (/ 1 (sqrt (- 1 (/ (* 2 phi) (* M M))))) (expt E phi))))

(define (round-off z n)
  (let ((power (expt 10 n)))
    (/ (round (* power z)) power)))

(define ACCURANCY 0.0000001)
(define MAXITERS 100)
;;calculates the root in a space [start, end] of f with deriv df
(define NewtonsMethod
	(lambda (f df start end)
		(letrec ((iter (lambda (g stats)
                                    (begin
                                      (format #t "iter g:~A fg:~A stat:~A~%" g (f g) stats)
		                      (if (< (* (- start g) (- g end)) 0)
 					'ERROUTBOUND
                                        (if (< (abs (nextx g)) ACCURANCY)
		      		          (list g stats (nextx g) )
				          (iter (nextx g) (+ stats 1)))))))
                         (nextx (lambda (x)
                                  (- x (/ (f x) (df x))))))
		  (format #t "starting~%")
                  (iter (* 0.5 (+ start end)) 0))))

(define (stats-iter-inc stats num)
  (+ stats num))

(define (stats-iter stats)
  stats)

;;calculates the root in a space [start, end] of f with deriv df
(define NewtonsSafeMethod
	(lambda (f df start end)
		(letrec ((dx-newton (lambda (x)
                           (/ (f x) (df x))))
                         (dx-bisect (lambda (xh xl)
                           (* 0.5 (- xh xl))))
                         (nextx-newton (lambda (x)
                           (- x (dx-newton x))))
                         (nextx-bisect (lambda (xh xl)
                           (+ xl (dx-bisect xh xl)))) 
        		 (iter (lambda (g xh xl dxold stats)
                                    (begin
                        ;              (format #t "iter g:~A fg:~A xh:~A xl:~A dxold:~A stat:~A~%" g (f g) xh xl dxold stats)
                                      (if (> (stats-iter stats) MAXITERS)
        			 	'ERRMAXITERS
                                         (let ((fg (f g))
                                               (dfg (df g)))
                                           (if (or (> (* (- (* (- g xh) dfg) fg) (- (* (- g xl) dfg) fg)) 0) ;;bisect if newton out of range
        					   (> (abs (* 2 fg)) (abs (* dxold dfg)))) ;;not decreasing fast enough
        				     (if (or (= xl (nextx-bisect xh xl)) (< (abs (nextx-bisect xh xl)) ACCURANCY)) ;; do bisect
                                               xl ;; no change or converge
                                               (let* ((g2 (nextx-bisect xh xl))
        					      (fg2 (f g2))
        					      (dfg2 (df g2)))
        					(if (< fg2 0)
        					   (iter g2 xh g2 (dx-bisect xh xl) (stats-iter-inc stats 1))
        					   (iter g2 g2 xl (dx-bisect xh xl) (stats-iter-inc stats 1)))))
                                             (if (or (= g (nextx-newton g)) (< (abs (nextx-newton g)) ACCURANCY)) ;; do newton
                                               g;; no change or converge
                                               (let* ((g2 (nextx-newton g))
        					      (fg2 (f g2))
        					      (dfg2 (df g2)))
        					(if (< fg2 0)
        					   (iter g2 xh g2 (dx-newton g) (stats-iter-inc stats 1))
        					   (iter g2 g2 xl (dx-newton g) (stats-iter-inc stats 1)))))))))))
        		(fl (f start))
        		(fh (f end))) 
		  (if (or (and (> fl 0) (> fh 0)) (and (< fl 0) (< fh 0)))
		    'ERRROOTMUSTBEBRACKETED) 
                  (if (= fl 0) start)
                  (if (= fh 0) end)
                  (if (< fl 0)
		    (iter (* 0.5 (+ start end)) end start (abs (- end start)) 0)
		    (iter (* 0.5 (+ start end)) start end (abs (- end start)) 0)))))
		  

(define SagdeevNewtonAtM
	(lambda (M start end)
                (let ((sAtM (lambda (phi) (sagdeev phi M)))
                      (dsAtM (lambda (phi) (DsagdeevDphi phi M))))
                  (NewtonsMethod sAtM dsAtM start end))))
             
 (define SagdeevNewtonSafeAtM
	(lambda (M start end)
                (let ((sAtM (lambda (phi) (sagdeev phi M)))
                      (dsAtM (lambda (phi) (DsagdeevDphi phi M))))
                  (NewtonsSafeMethod sAtM dsAtM start end))))

(define SagdeevRootsInSpace
	(lambda (start end step rA rB)
		(let* ((Ms (genlist start end step))
		       (sgdv (lambda (x) (SagdeevNewtonSafeAtM x rA rB)))
		       (roots (map sgdv Ms)))
                  (map cons Ms roots))))

(define variateSagdeevRoots
	(lambda (start end step rA rB rstep)
		(begin
		(format #t "foo~%")
		(letrec((rootvals (make-hash-table))
		      (iter (lambda (irA irB ht)
		              (if (< irB rB)
                                (begin (format #t "start:~A end:~A step:~A irA:~A irb: ~A rb:~A ~%" start end step irA irB rB)
                                (let ((vals (SagdeevRootsInSpace start end step irA irB)))
                                  (do-list i vals
                                    (if (and (not (eq? (cdr i) 'ERRMAXITERS)) (= (cdr i) irB))
				      (iter (+ irA rstep) (+ irB rstep) ht)
                                      (if (and (not (hash-table-exists? ht (car i))) (not (eq? (cdr i) 'ERRMAXITERS)))
                                        (hash-table-set! ht (car i) (cdr i)))))))
                                (begin (format #t "bar") 'ENDITER)))))
                  (iter rA (- rB 0.5) rootvals)
                  rootvals))))
                                  


;;(define (number-of-digits-gen x acc mod func)
;;  (if (= mod 0)
;;      acc
;;      (number-of-digits-gen (func x 10) (+ acc 1) (remainder (floor (func x 10)) 10) func)))

;;(define (number-of-digits x)
;;  (cond ((= x 1) 1)
;;        ((< x 1) (number-of-digits-gen x 0 1 *))
;;        ((> x 1) (number-of-digits-gen x 0 1 /))))

(define (number-of-digits x)
  (let ((parts (string-split (number->string x) ".")))
    (if (= 2 (length parts))
        (+ (string-length (car parts)) (string-length (cadr parts)))
        (string-length (car parts)))))

(define genlist
	(lambda (start end step)
		(if (> start end) '()
		    (cons start (genlist (round-off (+ start step) (number-of-digits step)) end step)))))

(define genpoints
        (lambda (start end amount)
                (letrec* ((dx (/ (- end start) amount))
                          (numdigits (number-of-digits dx))
                          (loop (lambda (s e) (if (> s e) '() (cons s (loop (round-off ( + s dx) numdigits) e))))))
                  (if (> start end)
                      '()
                      (loop start end)))))
                      
                    

(define linspace1 (genlist 0 1 0.0001))

(define res1 (map (lambda (s) (cons s (sagdeev s 1.5))) linspace1))

(define sagdeev4mach
	(lambda (M)
		(let* ((linspace (genlist 0 1 0.01))
		       (res (map (lambda (s) (cons s (sagdeev s M))) linspace)))
		  res)))

(define termsphiOfX
	(lambda (M phir dx terms)
          (letrec ((phiofx (lambda (i)
                           (cond ((hash-table-exists? ht (* i dx)) (hash-table-ref ht (* i dx))) 
				 ((= i 0) (let ((val  phir)) (format #t "1: ~A~%" val) (hash-table-set! ht (* i dx) val) val))
                                 ((= i 1) (let ((val (- phir (/ (* dx dx (DsagdeevDphi phir M)) 2)))) (format #t "-1: ~A~%" val) (hash-table-set! ht (* i dx) val) val))
                                 (else (let ((val (- (* 2 (phiofx (- i 1))) (phiofx (- i 2)) (* dx dx (DsagdeevDphi (phiofx (- i 1)) M))))) (format #t "~A: ~A~%" i val) (hash-table-set! ht (* i dx) val) val)))))
                   (ht (make-hash-table)))
            (phiofx terms)
            ht)))
          
(define symetricht 
	(lambda (ht)
	  (hash-table-for-each ht (lambda (k v) (hash-table-set! ht (- k) v)))))

(define write-plot
	(lambda (fname data)
		(with-output-to-file fname
			(lambda () (map (lambda (e) (format #t "~A    ~A~%" (car e) (cdr e))) data)))))

(define writeloop
	(lambda (M1 M2 step fname)
	  (let ((mspace (genlist M1 M2 step)))
            (do-list i mspace 
              (write-plot (string-append fname "-" (number->string i))
                          (sagdeev4mach i))))))

(define readrootsfile
	(lambda (fname)
	  (with-input-from-file fname
	    (lambda ()
              (letrec ((loop (lambda (lst)
			       (let ((line (read-line)))
                                 (if (eof-object? line)
                                   lst
                                   (let* ((strlist (string-split line " "))
                                          (numlist (map string->number strlist)))
                                     (loop (cons numlist lst))))))))
	        (loop '()))))))

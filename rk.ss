;;Stream computational approach for diff eq solving.
;;<dsp@2f30.org>
(declare (unit rk4))
;(module rk (integrate-system wave-eq do-n)
(import
  (rename scheme (force r5rs:force) (delay r5rs:delay))
  (rename chicken (promise? r5rs:promise?))
  )
(require-extension srfi-45)
;	 (require-extension srfi-41)
(use streams)

(define rk4
  (lambda  (f h)
    (let ((*h (scale-vector h))
	  (*2 (scale-vector 2))
	  (*1/2 (scale-vector (/ 1 2)))
	  (*1/6 (scale-vector (/ 1 6))))
      (lambda (y)
	(let*  ((k0 (*h (f y)))
		(k1 (*h (f (add-vectors y (*1/2 k0)))))
		(k2 (*h (f (add-vectors y (*1/2 k1)))))
		(k3 (*h (f (add-vectors y k2)))))
	  (add-vectors y
		       (*1/6 (add-vectors k0
					  (*2 k1)
					  (*2 k2)
					  k3))))))))

(define euler-1step
  (lambda (f h)
    (let ((*h (scale-vector h)))
      (lambda (y)
	(let* ((next (*h (f y))))
	  (add-vectors y next)
	  (vector-set! next 0 1))))))

(define integrate-system
  (lambda (sysderiv initstate h)
    ;		(let ((next (rk4 sysderiv h)))
    (let ((next (euler-1step sysderiv h)))
      (letrec ((states
		 (cons initstate
		       (delay (map-streams next
					   states)))))
	states))))

(define map-streams
  (lambda (f s)
    (cons (f (head s))
	  (delay (map-streams f (tail s))))))

(define head car)
(define tail
  (lambda (stream) (force (cdr stream))))

(define elemwise
  (lambda (f)
    (lambda vectors
      (generate-vector
	(vector-length (car vectors))
	(lambda (i)
	  (apply f
		 (map (lambda (v) (vector-ref v i))
		      vectors)))))))

(define generate-vector
  (lambda (size proc)
    (let ((ans (make-vector size)))
      (letrec ((loop
		 (lambda (i)
		   (cond ((= i size) ans)
			 (else
			   (vector-set! ans i (proc i))
			   (loop (+ i 1)))))))
	(loop 0)))))

(define add-vectors
  (elemwise +))

(define scale-vector
  (lambda (s)
    (elemwise (lambda (x) (* x s)))))

(define square
  (lambda (x)
    (* x x)))

;;sample diffeqs

(define damped-oscillator
  (lambda (R L C)
    (lambda (state)
      (let ((Vc (vector-ref state 0))
	    (Il (vector-ref state 1)))
	(vector  (- 0 (+ (/ Vc (* R C)) (/ Il C)))
		 (/ Vc L))))))

(define wave-eq
  (lambda (c)
    (lambda (state)
      (let ((ut (vector-ref state 0))
	    (ux (vector-ref state 1)))
	(vector (* c  ux)
		(- 0 (/ ut c)))))))

(define sampleode1
  (lambda ()
    (lambda (state)
      (let ((t (vector-ref state 0))
	    (ut (vector-ref state 1))
	    (ut1 (vector-ref state 2))
	    (ut2 (vector-ref state 3)))
	(vector (+ t 0.01)
		ut1
		ut2
		(- (* (sin t)
		      (sqrt (+ 1
			       (square ut2)))
		      ut1)
		   (/ ut
		      (+ 1
			 (exp (- 0 t))))
		   (- (square t))))))))


(define (do-n proc num s)
  (letrec ((loop
	     (lambda (p n s)
	       (cond ((< n num )
		      (begin
			(p (head s))
			(loop  p (+ n 1) (force (cdr s)))))))))
    (loop proc 1 s)))
;)

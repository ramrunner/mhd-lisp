(declare (uses rk4))
;(define the-states
;	(integrate-system
;		(wave-eq 1.2)
;		'#(1 0.1)
;		0.001))
;(define the-states
;	  (integrate-system
;			    (damped-oscillator 10000 1000 .001)
;					    '#(1 0)
;							    .01))

;print2file is a function
;that takes as input a list
;of columns and a filename
;and prints the corresponding columbs to it.

(define pr2f
  (lambda (what)
    (let ((cols '(0 1))
	  (fname "testdata"))
      (with-output-to-file fname (lambda ()
				   (if (not (null? what))
				     (let ((res (what)))
				       (if (not (null? res))
					 res))))))))

(define the-states
  (integrate-system
    (sampleode1)
    '#(0 1 2 3)
    .01))

(with-output-to-file "testdata"
		     (lambda ()
		       (do-n (lambda (i) (print (vector-ref i 0) 
						"\t" 
						(vector-ref i 1))) 
			     10 
			     the-states)))

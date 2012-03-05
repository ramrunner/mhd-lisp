; coordinates code by DsP
; TODO: implement type abstraction
;     : implement sytax rules for more than 2 ops
; a n-vector is a sequence of n nums
; make a vector: (mvec 1 2 3)
(define (mvec . list) list)
; make a matrix: (mmat '(1 2 3) '(4 5 6) '(7 8 9))
(define (mmat . list) list)

; sample data:
(define i3 (mmat '(1 0 0) '(0 1 0) '(0 0 1)))
(define backwards3 (mmat '(0 0 1) '(0 1 0) '(1 0 0)))
(define nine (mmat '(1 2 3) '(4 5 6) '(7 8 9)))
(define get-x (mmat '(1 0 0) '(0 0 0) '(0 0 0)))
(define v1 (mvec 1 2 3))

; mapcar
(define (mapcar func list)
  (if (null? list) '()
      (cons (func (car list)) (mapcar func (cdr list)))))

;vector add element by element
(define (v+ v1 v2)
  (map (lambda (a b) (+ a b)) v1 v2))
;vector subtract element by element
(define (v- v1 v2)
  (map (lambda (a b) (- a b)) v1 v2))

; matrix transpose
(define (mtranspose mat)
  (define (add-first-column vector matrix)
    (if (null? matrix) (mapcar (lambda (l) (list l)) vector)
	(cons (cons (car vector) (car matrix))
	      (add-first-column (cdr vector) (cdr matrix)))))
  (if (null? mat) '()
      (add-first-column (car mat) (mtranspose (cdr mat))))  )

; dot product
(define (mdot v1 v2)
  (if (null? v1) 0
      (+ (* (car v1) (car v2)) (mdot (cdr v1) (cdr v2)))))

; matrix multiply:
(define (mm m1 m2)
  (mapcar (lambda (rowfirstmat) 
	    (mapcar (lambda (columnsecondmat)
		      (mdot rowfirstmat columnsecondmat))
		    (mtranspose m2))) m1))

; scalar multiply
(define (mscmul s d)
  (mapcar (lambda (i) (if (pair? i) (mscmul s i) (* s i))) d))

; magnitude
(define (mmag v)
  (sqrt (mdot v v)))

(define (mnorm v)
  (mscmul (/ (mmag v)) v))

; a matrix-vector multiply, treating the vector as a column vector.
; (mmv (mmat '(1 2 3) '(4 5 6) '(7 8 9)) (mvec 1 2 3))
(define (mmv mat vec)
  (if (null? mat) '()
      (cons (mdot (car mat) vec) (mmv (cdr mat) vec))))


; this will make a disk of radius rad and density defined in e
; e is the "radius" density and d the "angular" density for
; creating the covering . n is the number of the gen loop runs
(define (mdisk rad e d n)
  (let loop ((a (mvec 0 0)) (res '()))
    (if (< (mmag a) rad)
      (begin (loop (mvec 0 (+ e (cadr a))) (cons a res)))
      (begin
        (let ((a (reverse res))
              (i2 (mmat '(1 0) '(0 1)))
              (R (mmat '(0 1) '(-1 0))))
          (let loopi ((ni 0) (b '()))
            (if (< ni n)
              (loopi (+ ni 1)
                     (cons (map (lambda (a) (v+ (mmv i2 a) (mmv (mscmul d R) a)))
                                (if (null? b) a (car b)))
                           b))
              b)))))))


(define write-list-to-file
  (lambda (filename ls)
    (with-output-to-file filename (lambda () (display-list ls)))))

;that was hard
(define display-list
  (lambda (ls)
    (if (null? ls)
      (newline)
      (begin
        (if (not (list? (car ls)))
          (begin
            (display (car ls))
            (display "\t")
            (display (cadr ls))
            (newline))
          (begin (display-list (car ls)) (display-list (cdr ls))))))))

(define yeah (mdisk 1 0.1 0.1 80))
(write-list-to-file "kota.dat" yeah)
yeah

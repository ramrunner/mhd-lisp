; coordinates code by DsP
; TODO: implement type abstraction
;     : implement sytax rules for more than 2 ops
; a n-vector is a sequence of n nums
; make a vector: (mvec 1 2 3)

;;LIBS
;; we want fold goddamit!!
(use srfi-1)

(define (mvec . list) list)
; make a matrix: (mmat '(1 2 3) '(4 5 6) '(7 8 9))
(define (mmat . list) list)

; sample data:
(define i3 (mmat '(1 0 0) '(0 1 0) '(0 0 1)))
(define backwards3 (mmat '(0 0 1) '(0 1 0) '(1 0 0)))
(define nine (mmat '(1 2 3) '(4 5 6) '(7 8 9)))
(define get-x (mmat '(1 0 0) '(0 0 0) '(0 0 0)))
(define v1 (mvec 1 2 3))

(define (caddddr p)
	(car (cdr (cdr (cdr (cdr p))))))

(define (cadddddr p)
	(car (cdr (cdr (cdr (cdr p))))))


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
;(define (mdisk rad e d n)
;  (let loop ((a (mvec 0 0)) (res '()))
;    (if (< (mmag a) rad)
;      (begin (loop (mvec 0 (+ e (cadr a))) (cons a res)))
;      (begin
;        (let ((a (reverse res))
;              (i2 (mmat '(1 0) '(0 1)))
;              (R (mmat '(0 1) '(-1 0))))
;          (let loopi ((ni 0) (b '()))
;            (if (< ni n)
;              (loopi (+ ni 1)
;                     (cons (map (lambda (a) (v+ (mmv i2 a) (mmv (mscmul d R) a)))
;                                (if (null? b) a (car b)))
;                           b))
;              b)))))))
;; with rot matrix
(define (mdisk rad e d n)
  (let loop ((a (mvec 0 0)) (res '()))
    (if (< (mmag a) rad)
      (begin (loop (mvec 0 (+ e (cadr a))) (cons a res)))
      (begin
	(let ((a (reverse res))
	      (rot (lambda (ang) 
		     (let* ((x (eval (cos ang)))
			    (y (eval (sin ang)))
			    (z (eval (- y))))
		       (list (list x z) (list y x))))))
	  (let loopi ((ni 0) (b '()))
	    (if (< ni n)
	      (loopi (+ ni (/ 3.14159 180))
		     (cons (map (lambda (a) (mmv (rot ni) a))
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

;; XXX add types and checks
;; a line is created from vectors
(define get-line
  (lambda (ls)
    (if (not (null? ls))
      (car ls))))

;;vectors are created from points
(define get-vector
  (lambda (ls)
    (if (not (null? ls))
      (car ls))))

;(define print-radii
;  (lambda (ls)
;    (if (not (null? ls))
;      (begin
;	(let ((allvecs (map (lambda (ls) ls) ls))
;; smth like that?
;;(map (lambda (ls) (map (lambda (ls) (fold (lambda (a b) (+ a b)) (car ls) (cdr ls))) ls) yeah)

(define square (lambda (x) (* x x)))

(define display-list3d
  (lambda (ls)
    (if (null? ls)
      (newline)
      (begin
	(if (not (list? (car ls)))
	  (begin
	    (display (car ls))
	    (display "\t")
	    (display (cadr ls))
	    (display "\t")
	    (display (caddr ls))
	    (newline))
	  (begin (display-list3d (car ls)) (display-list3d (cdr ls))))))))

(define write-list3d-to-file
  (lambda (filename ls)
    (with-output-to-file filename (lambda () (display-list3d ls)))))

(define yeah (mdisk 1 0.1 0.1 3.14159))
;(write-list-to-file "kota.dat" yeah)
yeah
;(map (lambda (r) (map (lambda (p) (+ (square (car p))  (square (cadr p)))) r)) yeah)
(define points (map (lambda (r)
		      (map (lambda (p)
			     (let* ((radsq (+ (square (car p)) (square (cadr p))))
				    (x1 (/ (car p) (+ 1 (/ radsq 4))))
				    (x2 (/ (cadr p) (+ 1 (/ radsq 4))))
				    (x3 (/ (- 1 (/ radsq 4)) (+ 1 (/ radsq 4)))))
			       (list x1 x2 x3)))
			   r))
		    yeah))
;(write-list3d-to-file "3d.dat" points)

(define (p3d x y z #!optional (e i3) (b i3) (di (mdisk 1 0.1 1.3 3)))
  (list x y z e b di))

(define (p3d-get-pos p)
  (list (car p) (cadr p) (caddr p)))

(define (p3d-get-e p)
  (cadddr p))

(define (p3d-get-b p)
  (caddddr p))

(define (p3d-get-disk p)
  (car (cdr (cddddr p))))

(define plist (list (p3d 0 0 0) (p3d 0 1 0) (p3d 1 1 1)))

(define (print-p3d p)
  (display "3dPoint @ ")
  (display (p3d-get-pos p))
  (display " Efield: ")
  (display (p3d-get-e p))
  (display " Bfield: ")
  (display (p3d-get-b p))
  (newline)
  (display "Disk:")
  (display (p3d-get-disk p))
  (newline))

(map print-p3d plist)
(newline)

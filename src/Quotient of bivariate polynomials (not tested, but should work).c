#include "theta_implementation.h"

struct polynomial polysum(struct polynomial p1,struct polynomial p2){
    struct polynomial output;
    int counter;
    if (p1.degree>p2.degree){
        output.degree = p1.degree;
        output.coeffs = malloc((output.degree+1)*sizeof(float));
        for (counter = 0; counter <= p2.degree; counter++){
            output.coeffs[counter] = p1.coeffs[counter]+p2.coeffs[counter];
        }
        for (counter = p2.degree; counter <= p1.degree; counter++){
            output.coeffs[counter] = p1.coeffs[counter];
        }
    }
    else{
        output.degree = p2.degree;
        output.coeffs = malloc((output.degree+1)*sizeof(float));
        for (counter = 0; counter <= p1.degree; counter++){
            output.coeffs[counter] = p1.coeffs[counter]+p2.coeffs[counter];
        }
        for (counter = p1.degree; counter <= p2.degree; counter++){
            output.coeffs[counter] = p2.coeffs[counter];
        }
    }
    return output;
}
struct polynomial polydiff(struct polynomial p1,struct polynomial p2){
    struct polynomial output;
    int counter;
    if (p1.degree>p2.degree){
        output.degree = p1.degree;
        output.coeffs = malloc((output.degree+1)*sizeof(float));
        for (counter = 0; counter <= p2.degree; counter++){
            output.coeffs[counter] = p1.coeffs[counter]-p2.coeffs[counter];
        }
        for (counter = p2.degree; counter <= p1.degree; counter++){
            output.coeffs[counter] = p1.coeffs[counter];
        }
    }
    else if (p2.degree>p1.degree){
        output.degree = p2.degree;
        output.coeffs = malloc((output.degree+1)*sizeof(float));
        for (counter = 0; counter <= p1.degree; counter++){
            output.coeffs[counter] = p1.coeffs[counter]-p2.coeffs[counter];
        }
        for (counter = p1.degree; counter <= p2.degree; counter++){
            output.coeffs[counter] = -p2.coeffs[counter];
        }
    }
    else{
        output.degree = p1.degree;
        for(counter = p1.degree; counter>= 0; counter--){
            if (p1.coeffs[counter] == p2.coeffs[counter]){
                output.degree -= 1;
            }
            else{
                break;
            }
        }
        output.coeffs = malloc((output.degree+1)*sizeof(float));
        for (counter = 0; counter <= output.degree; counter++){
            output.coeffs[counter] = p1.coeffs[counter]-p2.coeffs[counter];
        }
    }
    return output;
}
struct polynomial polyprod(struct polynomial p1, struct polynomial p2){
    struct polynomial output;
    output.degree = p1.degree + p2.degree;
    output.coeffs = malloc((output.degree+1)*sizeof(float));
    int counter;
    int counter1;
    for (counter = 0; counter<= output.degree; counter++){
        output.coeffs[counter] = 0;
        for (counter1 = 0; counter1<= counter; counter++){
            if (counter1<= p1.degree && counter-counter1<=p2.degree){
                output.coeffs[counter] += p1.coeffs[counter1]*p2.coeffs[counter-counter1];
            }
        }
    }
    return output;
}
struct polynomial polydiv(struct polynomial p1, struct polynomial p2){//does not give remainder.
    struct polynomial output;
    struct polynomial p1mod;
    if(p1.degree<p2.degree){
        output.degree = 0;
        output.coeffs = malloc(sizeof(float));
        output.coeffs[0] = 0;
        return output;
    }
    else{
        output.degree = p1.degree-p2.degree;
        output.coeffs = malloc((output.degree+1)*sizeof(float));
        int counter;
        for (counter = 0; counter< output.degree; counter++){
            output.coeffs[counter] = 0;
        }
        output.coeffs[output.degree] = p1.coeffs[p1.degree]/p2.coeffs[p2.degree];
        return polysum(output,polydiv(polydiff(p1,polyprod(output,p2)),p2));
    }
}
int eval(struct polynomial p1, float n){
    int output = 0;
    int counter;
    for (counter = 0; counter<= p1.degree; counter++){
        output = n*output + p1.coeffs[counter];
    }
    return output;
}
struct rational{
struct polynomial numerator;
struct polynomial denominator;
};
struct rational ratsum(struct rational r1, struct rational r2){
    struct rational output;
    output.numerator.degree = polysum(polyprod(r1.numerator, r2.denominator), polyprod(r1.denominator, r2.numerator)).degree;
    output.numerator.coeffs = malloc((output.numerator.degree+1)*sizeof(float));
    int counter;
    for (counter = 0; counter<=output.numerator.degree; counter++){
        output.numerator.coeffs[counter] = polysum(polyprod(r1.numerator, r2.denominator), polyprod(r1.denominator, r2.numerator)).coeffs[counter];
    }
    output.numerator.degree = polyprod(r1.denominator, r2.denominator).degree;
    output.numerator.coeffs = malloc((output.numerator.degree+1)*sizeof(float));
    for (counter = 0; counter<=output.numerator.degree; counter++){
        output.numerator.coeffs[counter] = polyprod(r1.denominator, r2.denominator).coeffs[counter];
    }
    return output;
}
struct rational ratdiff(struct rational r1, struct rational r2){
    struct rational output;
    output.numerator.degree = polydiff(polyprod(r1.numerator, r2.denominator), polyprod(r1.denominator, r2.numerator)).degree;
    output.numerator.coeffs = malloc((output.numerator.degree+1)*sizeof(float));
    int counter;
    for (counter = 0; counter<=output.numerator.degree; counter++){
        output.numerator.coeffs[counter] = polydiff(polyprod(r1.numerator, r2.denominator), polyprod(r1.denominator, r2.numerator)).coeffs[counter];
    }
    output.numerator.degree = polyprod(r1.denominator, r2.denominator).degree;
    output.numerator.coeffs = malloc((output.numerator.degree+1)*sizeof(float));
    for (counter = 0; counter<=output.numerator.degree; counter++){
        output.numerator.coeffs[counter] = polyprod(r1.denominator, r2.denominator).coeffs[counter];
    }
    return output;
}
struct rational ratprod(struct rational r1, struct rational r2){
    struct rational output;
    output.numerator.degree = polyprod(r1.numerator, r2.numerator).degree;
    output.numerator.coeffs = malloc((output.numerator.degree+1)*sizeof(float));
    int counter;
    for (counter = 0; counter<=output.numerator.degree; counter++){
        output.numerator.coeffs[counter] = polyprod(r1.numerator, r2.numerator).coeffs[counter];
    }
    output.numerator.degree = polyprod(r1.denominator, r2.denominator).degree;
    output.numerator.coeffs = malloc((output.numerator.degree+1)*sizeof(float));
    for (counter = 0; counter<=output.numerator.degree; counter++){
        output.numerator.coeffs[counter] = polyprod(r1.denominator, r2.denominator).coeffs[counter];
    }
    return output;
}
struct rational ratquot(struct rational r1, struct rational r2){
    struct rational output;
    output.numerator.degree = polyprod(r1.numerator, r2.denominator).degree;
    output.numerator.coeffs = malloc((output.numerator.degree+1)*sizeof(float));
    int counter;
    for (counter = 0; counter<=output.numerator.degree; counter++){
        output.numerator.coeffs[counter] = polyprod(r1.numerator, r2.denominator).coeffs[counter];
    }
    output.numerator.degree = polyprod(r1.denominator, r2.numerator).degree;
    output.numerator.coeffs = malloc((output.numerator.degree+1)*sizeof(float));
    for (counter = 0; counter<=output.numerator.degree; counter++){
        output.numerator.coeffs[counter] = polyprod(r1.denominator, r2.numerator).coeffs[counter];
    }
    return output;
}
struct squaremat{
int size; //degree in each
float **coeffs;
};
struct tvalues{
int size;
float *values;
};
struct squaremat bivarprod(struct polynomial pt1, struct polynomial pt2){//we assume same degree b/c the thing is symmetric in T1,T2
    struct squaremat output;
    output.size = pt1.degree;
    output.coeffs = malloc((output.size+1)*sizeof(float*));
    int counter;
    int counter1;
    for (counter = 0; counter<= output.size; counter++){
        output.coeffs[counter] = malloc((output.size+1)*sizeof(float));
    }
    for (counter = 0; counter<= output.size; counter++){
        for (counter1 = 0; counter1<= output.size; counter1++){
            output.coeffs[counter][counter1] = pt1.coeffs[counter]*pt2.coeffs[counter1];
        }
    }
    return output;
}
struct polynomial bigprod(struct tvalues tvals){
    if (tvals.size == 0){
        struct polynomial output;
        output.degree = 0;
        output.coeffs = malloc(sizeof(float));
        output.coeffs[0] = 1;
        return output;
    }
    else{
        struct tvalues tmodvals;
        tmodvals.size = tvals.size-1;
        tmodvals.values = malloc((tmodvals.size+1)*sizeof(float));
        int counter;
        for (counter = 0; counter<= tmodvals.size; counter++){
            tmodvals.values[counter] = tvals.values[counter];
        }
        struct polynomial last;
        last.degree = 1;
        last.coeffs = malloc(2*sizeof(float));
        last.coeffs[0] = -tvals.values[tvals.size];
        last.coeffs[1] = 1;
        return polyprod(bigprod(tmodvals),last);
    }
}
struct squaremat squaresum(struct squaremat s1, struct squaremat s2){//assumes same size. Can modify
    struct squaremat output;
    output.size = s1.size;
    output.coeffs = malloc((output.size+1)*sizeof(float*));
    int counter;
    int counter1;
    for (counter = 0; counter<= output.size; counter++){
        output.coeffs[counter] = malloc((output.size+1)*sizeof(float));
    }
    for (counter = 0; counter<= output.size; counter++){
        for (counter1 = 0; counter1<= output.size; counter1++){
            output.coeffs[counter][counter1] = s1.coeffs[counter][counter1]+s2.coeffs[counter][counter1];
        }
    }
    return output;
}
struct squaremat squarediff(struct squaremat s1, struct squaremat s2){//Assumes same size. Can modify.
    struct squaremat output;
    output.size = s1.size;
    output.coeffs = malloc((output.size+1)*sizeof(float*));
    int counter;
    int counter1;
    int counter2 = 0;
    for (counter = 0; counter<= output.size; counter++){
        output.coeffs[counter] = malloc((output.size+1)*sizeof(float));
    }
    for (counter = 0; counter<= output.size; counter++){
        for (counter1 = 0; counter1<= output.size; counter1++){
            output.coeffs[counter][counter1] = s1.coeffs[counter][counter1]-s2.coeffs[counter][counter1];
            if (counter == output.size || counter1 == output.size){
                if (output.coeffs[counter][counter1] == 0){
                    counter2++;
                }
            }
        }
    }
    if (counter2++ == 2*output.size + 1){
        struct squaremat reoutput;
        reoutput.size = output.size-1;
        reoutput.coeffs = malloc((reoutput.size+1)*sizeof(float*));
        int counter;
        int counter1;
        int counter2 = 0;
        for (counter = 0; counter<= reoutput.size; counter++){
            output.coeffs[counter] = malloc((reoutput.size+1)*sizeof(float));
        }
        for (counter = 0; counter<= reoutput.size; counter++){
            for (counter1 = 0; counter1<= reoutput.size; counter1++){
                reoutput.coeffs[counter][counter1] = output.coeffs[counter][counter1];
            }
        }
        return reoutput;
    }
    else{
    return output;
    }
}
struct squaremat squareprod(struct squaremat s1, struct squaremat s2){//Assumes same size. Can modify.
    struct squaremat output;
    output.size = s1.size+s2.size;
    output.coeffs = malloc((output.size+1)*sizeof(float*));
    int counter;
    int counter1;
    int counter2;
    int counter3;
    for (counter = 0; counter<= output.size; counter++){
        output.coeffs[counter] = malloc((output.size+1)*sizeof(float));
    }
    for (counter = 0; counter<= output.size; counter++){
        for (counter1 = 0; counter1<= output.size; counter1++){
            output.coeffs[counter][counter1] = 0;
        }
    }
    for (counter = 0; counter<= s1.size; counter++){
        for (counter1 = 0; counter1<= s1.size; counter1++){
            for (counter2 = 0; counter2<= s2.size; counter2++){
                for (counter3 = 0; counter3<= s2.size; counter3++){
                    output.coeffs[counter+counter2][counter1+counter3] = output.coeffs[counter+counter2][counter1+counter3]+ s1.coeffs[counter][counter1]*s2.coeffs[counter2][counter3];    
                }
            }
        }
    }
    return output;
}
struct squaremat scalardiv(struct squaremat s1, float n){
    struct squaremat output;
    output.size = s1.size;
    output.coeffs = malloc((output.size+1)*sizeof(float*));
    int counter;
    int counter1;
    for (counter = 0; counter<= output.size; counter++){
        output.coeffs[counter] = malloc((output.size+1)*sizeof(float));
    }
    for (counter = 0; counter<= output.size; counter++){
        for (counter1 = 0; counter1<= output.size; counter1++){
            output.coeffs[counter][counter1] = s1.coeffs[counter][counter1]/n;
        }
    }
    return output;
}
struct squaremat lagrange(struct tvalues t1vals, struct tvalues t2vals, struct squaremat funcvals){//t1 vals and t2 vals are the values I put in for t1, t2 in f(t1,t2), and funcvals is an array that gives the output.
    struct squaremat output;
    struct polynomial tempol1;//will be changed (temporary polynomial)
    tempol1.degree = 1;
    tempol1.coeffs = malloc(2*sizeof(float));
    tempol1.coeffs[1] = 1;
    struct polynomial tempol2;//will be changed (temporary polynomial)
    tempol2.degree = 1;
    tempol2.coeffs = malloc(2*sizeof(float));
    tempol2.coeffs[1] = 1;
    output.size = funcvals.size;
    output.coeffs = malloc((output.size+1)*sizeof(float*));
    int counter;
    int counter1;
    for (counter = 0; counter<= output.size; counter++){
        output.coeffs[counter] = malloc((output.size+1)*sizeof(float));
    }
    for (counter = 0; counter<= output.size; counter++){
        for (counter1 = 0; counter1<= output.size; counter1++){
            output.coeffs[counter][counter1] = 0;
        }
    }
    for (counter = 0; counter<= output.size; counter++){
        for (counter1 = 0; counter1<= output.size; counter1++){
            tempol1.coeffs[0] = -t1vals.values[counter];
            tempol2.coeffs[0] = -t2vals.values[counter1];
            output = squaresum(output,scalardiv(bivarprod(polydiv(bigprod(t1vals),tempol1),polydiv(bigprod(t2vals),tempol2)), eval(polydiv(bigprod(t1vals),tempol1),t1vals.values[counter])*eval(polydiv(bigprod(t2vals),tempol2),t2vals.values[counter1])/funcvals.coeffs[counter][counter1]));
        }
    }
    return output;
}
struct squareratmat{
int size;
struct rational **entries;
};
struct squaremat squarequot(struct squaremat s1, struct squaremat s2){
    struct squaremat temp;//Temporary; recursive algorithm
    if (s1.size<s2.size){
        temp.size = 1;
        temp.coeffs = malloc(sizeof(float*));
        temp.coeffs[0] = malloc(sizeof(float));
        return temp;
    }
    temp.size = s1.size-s2.size;
    temp.coeffs = malloc((temp.size+1)*sizeof(float*));
    int counter;
    int counter1;
    for (counter = 0; counter<= temp.size; counter++){
        temp.coeffs[counter] = malloc((temp.size+1)*sizeof(float));
    }
    int biggest1 = 0;
    int biggest2 = 0;
    for (counter = 0; counter<= s1.size; counter++){
        for (counter1 = 0; counter1<= s1.size; counter1++){
            if (s1.coeffs[counter][counter1]!= 0 && counter+counter1 > biggest1 + biggest2){
                biggest1 = counter;
                biggest2 = counter1;
            }    
        }
    }
    int biggest3 = 0;
    int biggest4 = 0;
    for (counter = 0; counter<= s2.size; counter++){
        for (counter1 = 0; counter1<= s2.size; counter1++){
            if (s1.coeffs[counter][counter1]!= 0 && counter+counter1 > biggest3 + biggest4){
                biggest3 = counter;
                biggest4 = counter1;
            }    
        }
    }
    for (counter = 0; counter<= temp.size; counter++){
        for (counter1 = 0; counter1<= temp.size; counter1++){
             temp.coeffs[counter][counter1]= 0;
        }
    }
    temp.coeffs[biggest1-biggest3][biggest2-biggest4] = s1.coeffs[biggest1][biggest2]/s2.coeffs[biggest3][biggest4];
    return squaresum(temp, squarediff(s1, squareprod(s2,temp)));
}

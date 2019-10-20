SAMPL ING D ISTR IBUT IONS
AND CONF IDENCE
In Chapters 15 and 16, you applied
the idea of a probability distribution
to examples where a random variable is
defined as some measurement or observation
of interest. In this chapter, you'll consider sample
statistics themselves as random variables to introduce
the concept of a sampling distribution-a probability distribution that is used
to account for the variability naturally present when you estimate population parameters using sample statistics. I'll then introduce the idea of a
confidence interval, which is a direct reflection of the variability in a sampling
distribution, used in a way that results in an interval estimate of a population
parameter. This will form the foundation for formal hypothesis testing in
Chapter 18.
17.1 Sampling Distributions
A sampling distribution is just like any other probability distribution, but
it is specifically associated with a random variable that is a sample statistic.
In Chapters 15 and 16, we assumed we knew the parameters of the relevant
example distribution (for example, the mean and the standard deviation
                      of a normal distribution or the probability of success in a binomial distribution), but in practice these kinds of quantities are often unknown. In these
cases, you'd typically estimate the quantities from a sample (see Figure 13-2
                                                              on page 266 for a visual illustration of this). Any statistic estimated from a
sample can be treated as a random variable, with the estimated value itself as
the realization of that random variable. It's therefore entirely possible that
different samples from the same population will provide a different value for
the same statistic-realizations of random variables are naturally subject to
variability. Being able to understand and model this natural variability inherent in estimated sample statistics (using relevant sampling distributions) is a
key part of many statistical analyses.
Like any other probability distribution, the central "balance" point of a
sampling distribution is its mean, but the standard deviation of a sampling
distribution is referred to as a standard error. The slight change in terminology reflects the fact that the probabilities of interest are no longer tied to
raw measurements or observations per se, but rather to a quantity calculated
from a sample of such observations. The theoretical formulas for various
sampling distributions therefore depend upon (a) the original probability
distributions that are assumed to have generated the raw data and (b) the
size of the sample itself.
This section will explain the key ideas and provide some examples, and
I'll focus on two simple and easily recognized statistics: a single sample mean
and a single sample proportion. I'll then expand on this in Chapter 18 when
covering hypothesis testing, and you'll need to understand the role of sampling distributions in assessing important model parameters when you look
at regression methods in Chapters 20 to 22.
NOTE The validity of the theory of sampling distributions as discussed in this chapter makes
an important assumption. Whenever I talk about a sample of data from which a given
statistic is calculated, I assume those observations are independent of one another and
that they are identically distributed. You'll see this notion-independent, identically
distributed observations-frequently abbreviated as iid in statistical material.
17.1.1 Distribution for a Sample Mean
The arithmetic mean is arguably the most common measure of centrality
(Section 13.2.1) used when summarizing a data set.
Mathematically, the variability inherent in an estimated sample mean
is described as follows: Formally, denote the random variable of interest as
X¯ . This represents the mean of a sample of n observations from the "raw
observation" random variable X, as in x1, x2,. . . , xn. Those observations are
assumed to have a true finite mean ?????? < µX < ??? and a true finite standard
deviation 0 < ??X < ???. The conditions for finding the probability distribution of a sample mean vary depending on whether you know the value of the
standard deviation.
368 Chapter 17
Situation 1: Standard Deviation Known
When the true value of the standard deviation ??X is known, then the following are true:
  . If X itself is normal, the sampling distribution of X¯ is a normal distribution, with mean µX and standard error ??X/
  ???
n.
. If X is not normal, the sampling distribution of X¯ is still approximately
normal, with mean µX and standard error ??X/
  ???
n, and this approximation improves arbitrarily as n ??? ???. This is known as the central limit
theorem (CLT).
Situation 2: Standard Deviation Unknown
In practice, you commonly won't know the true value of the standard deviation of the raw measurement distribution that generated your sample data.
In this eventuality, it's usual to just replace ??X with sX, which is the standard deviation of the sampled data. However, this substitution introduces
additional variability that affects the distribution associated with the sample
mean random variable.
. Standardized values (Section 16.2.2) of the sampling distribution of X¯
follow a t-distribution with ?? = n ??? 1 degrees of freedom; standardization
is performed using the standard error sX/
  ???
n.
. If, additionally, n is small, then it is necessary to assume the distribution
of X is normal for the validity of this t-based sampling distribution of X¯ .
The nature of the sampling distribution of X¯ therefore depends upon
whether the true standard deviation of the observations is known, as well
as the sample size n. The CLT states that normality occurs even if the raw
observation distribution is itself not normal, but this approximation is less
reliable if n is small. It's a common rule of thumb to rely on the CLT only
if n ??? 30. If sX, the sample standard deviation, is used to calculate the
standard error of X¯ , then the sampling distribution is the t-distribution
(following standardization). Again, this is generally taken to be reliable
only if n ??? 30.
Example: Dunedin Temperatures
As an example, suppose that the daily maximum temperature in the month
of January in Dunedin, New Zealand, follows a normal distribution, with a
mean of 22 degrees Celsius and a standard deviation of 1.5 degrees. Then,
in line with the comments for situation 1, for samples of size n = 5, the sampling distribution of X¯ will be normal, with mean 22 and standard error
1.5/
  ???
5 ??? 0.671.
The top image of Figure 17-1 shows the raw measurement distribution
along with this sampling distribution. You can produce this with code that's
familiar from Chapter 16.
Sampling Distributions and Confidence 369
R> xvals <- seq(16,28,by=0.1)
R> fx.samp <- dnorm(xvals,22,1.5/sqrt(5))
R> plot(xvals,fx.samp,type="l",lty=2,lwd=2,xlab="",ylab="")
R> abline(h=0,col="gray")
R> fx <- dnorm(xvals,22,1.5)
R> lines(xvals,fx,lwd=2)
R> legend("topright",legend=c("raw obs. distbn.","sampling distbn. (mean)"),
          lty=1:2,lwd=c(2,2),bty="n")
In this example, the sampling distribution of X¯ is clearly a taller, skinnier normal distribution than the one tied to the observations. This makes
sense-you expect less variability in an average of several measurements as
opposed to the raw, individual measurements. Furthermore, the presence
of n in the denominator of the standard error dictates a more precise distribution around the mean if you increase the sample size. Again, this makes
sense-means will "vary less" between samples of a larger size.
You can now ask various probability questions; note that distinguishing
between the measurement distribution and the sampling distribution is
important. For example, the following code provides Pr(X < 21.5), the probability that a randomly chosen day in January has a maximum temperature
of less than 21.5 degrees:
  R> pnorm(21.5,mean=22,sd=1.5)
[1] 0.3694413
The next bit of code provides the probability that the sample mean will
be less than 21.5 degrees, Pr(X¯ < 21.5), based on a sample of five random
days in January:
  R> pnorm(21.5,mean=22,sd=1.5/sqrt(5))
[1] 0.2280283
The line-shaded areas on the top of Figure 17-1 show these two probabilities. In R, these shaded areas can be added to that plot by running the
following lines directly after the earlier code:
  R> abline(v=21.5,col="gray")
R> xvals.sub <- xvals[xvals<=21.5]
R> fx.sub <- fx[xvals<=21.5]
R> fx.samp.sub <- fx.samp[xvals<=21.5]
R> polygon(cbind(c(21.5,xvals.sub),c(0,fx.sub)),density=10)
R> polygon(cbind(c(21.5,xvals.sub),c(0,fx.samp.sub)),density=10,
           angle=120,lty=2)
Note that in previous uses of polygon, you've simply specified a col; in this
example, I implemented shading lines instead, using the arguments density
(number of lines per inch) and angle (slope of lines in degrees; defaults to
                                      angle=45).
370 Chapter 17
Figure 17-1: Illustrating the sampling distribution of a sample mean
for n = 5, based on an N(22,1.5) raw observation distribution. Top:
  the normal-based version of the sampling distribution (assuming
                                                         ??X is known) compared to the observation distribution. Bottom:
  the t-based version of the sampling distribution, using 4 degrees
of freedom (in other words, assuming s has been used to calculate
            the standard error), compared to a standard normal. Shaded areas
represent Pr(X < 21.5), Pr(X¯ < 21.5) (solid and dashed, topmost
                                       plot) and Pr(T < (21.5 ??? x¯)/(s/
                                                                       ???
                                                                     5)) (dotted, bottom plot).
To evaluate the probabilities, note that you've required knowledge of
the parameters governing X. In practice, you'll rarely have these quantities
(as noted in situation 2). Instead, you obtain a sample of data and calculate
summary statistics.
Sampling Distributions and Confidence 371
Running the following line produces five randomly generated Dunedin
temperatures from the X ??? N(22,1.5) distribution:
  R> obs <- rnorm(5,mean=22,sd=1.5)
R> obs
[1] 22.92233 23.09505 20.98653 20.10941 22.33888
Now, for the sake of the example, say these five values constitute all the
data you have for this particular problem; in other words, pretend you don't
know that µX = 22 and ??X = 1.5. Your best guesses of the true values of µX
and ??X, denoted x¯ and s, respectively, are therefore as follows:
  R> obs.mean <- mean(obs)
R> obs.mean
[1] 21.89044
R> obs.sd <- sd(obs)
R> obs.sd
[1] 1.294806
The estimated standard error can be calculated:
  R> obs.mean.se <- obs.sd/sqrt(5)
R> obs.mean.se
[1] 0.5790549
Because n = 5 is relatively small, you must assume the values in obs are
realizations from a normal distribution, in line with the points made for situation 2. This allows you to handle the sampling distribution of X¯ using the
t-distribution with 4 degrees of freedom. Recall from Section 16.2.3, though,
that any t-distribution is typically placed on a standardized scale. Therefore,
for you to find the probability that the mean temperature (in a sample of
                                                           five days) is less than 21.5 based on your calculated sample statistics, you
must first standardize this value using the rules outlined in Section 16.2.2.
Label the corresponding random variable as T and the specific value as t4,
stored as the object t4 in R.
R> t4 <- (21.5-obs.mean)/obs.mean.se
R> t4
[1] -0.6742706
This has placed the value of interest, 21.5, on the standardized scale,
making it interpretable with respect to a standard normal distribution or,
as is correct in this setting (because you are using the estimate s rather than
                               the unknown ??X in calculating the standard error), t4 follows the aforementioned t-distribution with 4 degrees of freedom. The estimated probability is
as follows.
372 Chapter 17
R> pt(t4,df=4)
[1] 0.26855
Note that when you calculated the "true" theoretical probability from
the sampling distribution of Pr(X¯ < 21.5), you got a result of about 0.23 (see
                                                                            page 370), but the same probability based on standardization using sample
statistics of the data obs (in other words, estimates of the true theoretical
                            values Pr(T < t4)) has been computed as 0.27 (2 d.p.).
The bottom image of Figure 17-1 provides the t-distribution with ?? = 4,
marking off the probability described. The N(0,1) density is also plotted for
comparison; this represents the standardized version of the N(22,1.5/
                                                                ???
                                                              5)
sampling distribution from earlier, in situation 1. You can produce this
image with the following lines:
  R> xvals <- seq(-5,5,length=100)
R> fx.samp.t <- dt(xvals,df=4)
R> plot(xvals,dnorm(xvals),type="l",lty=2,lwd=2,col="gray",xlim=c(-4,4),
        xlab="",ylab="")
R> abline(h=0,col="gray")
R> lines(xvals,fx.samp.t,lty=3,lwd=2)
R> polygon(cbind(c(t4,-5,xvals[xvals<=t4]),c(0,0,fx.samp.t[xvals<=t4])),
           density=10,lty=3)
R> legend("topright",legend=c("N(0,1) standard","t (4 df)"),
          col=c("gray","black"),lty=2:3,lwd=c(2,2),bty="n")
Consideration of probability distributions associated with sample means
is clearly not a trivial exercise. Using sample statistics governs the nature
of the sampling distribution; in particular, it will be t based if you use the
sample standard deviation to calculate the standard error. However, as the
examples here have shown, once that's been established, the calculation of
various probabilities is easy and follows the same general rules and R functionality detailed in Section 16.2.
17.1.2 Distribution for a Sample Proportion
Sampling distributions for sample proportions are interpreted in much the
same way. If n trials of a success/failure event are performed, you can obtain
an estimate of the proportion of successes; if another n trials are performed,
the new estimate could vary. It's this variability that you're investigating.
The random variable of interest, P^, represents the estimated proportions of successes over any n trials, each resulting in some defined binary
outcome. It is estimated as p^ =
  x
n
, where x is the number of successes in a
sample of size n. Let the corresponding true proportion of successes (often
                                                                      unknown) simply be denoted with ??.
NOTE Note that ?? as used in this setting doesn't refer to the common geometric value 3.14
(2 d.p.). Rather, it's simply standard notation to refer to a true population proportion
using the ?? symbol.
Sampling Distributions and Confidence 373
The sampling distribution of P^ is approximately normal with mean ??
and standard error ???
??(1 ??? ??)/n. The following are the key things to note:
  . This approximation is valid if n is large and/or ?? is not too close to
either 0 or 1.
. There are rules of thumb to determine this validity; one such rule is to
assume the normal approximation is satisfactory if both n?? and n(1 ??? ??)
are greater than 5.
. When the true ?? is unknown or is unassumed to be a certain value, it is
typically replaced by p^ in all of the previous formulas.
As long as you can deem the approximation to the normal distribution
valid, this is the only probability distribution that you need to be concerned
with. However, it's worth noting that the standard error of the sampling distribution of a sample proportion depends directly upon the proportion ??.
This becomes important when constructing confidence intervals and carrying out hypothesis tests, which you'll begin to explore in Chapter 18.
Let's look at a practical example. Suppose a political commentator in
the United States is interested in the proportion of voting-age citizens in
her home city that already know how they will vote in the next presidential
election. She obtains a yes or no answer from 118 suitable randomly selected
individuals. Of these individuals, 80 say they know how they'll vote. To investigate the variability associated with the proportion of interest, you'll therefore need to consider
P^ ??? N *
  ,
p^,
r
p^(1 ??? p^)
n
+
  -
  , (17.1)
where p^ =
  80
118 . In R, the following gives you the estimate of interest:
  R> p.hat <- 80/118
R> p.hat
[1] 0.6779661
In the sample, about 68 percent of the surveyed individuals know how
they will vote in the next election. Note also that, according to the aforementioned rule of thumb, the approximation to the normal distribution is
valid because both values are greater than 5.
R> 118*p.hat
[1] 80
R> 118*(1-p.hat)
[1] 38
Estimate the standard error with the following:
  R> p.se <- sqrt(p.hat*(1-p.hat)/118)
R> p.se
[1] 0.04301439
374 Chapter 17
Then, you can plot the corresponding sampling distribution using
this code:
  R> pvals <- seq(p.hat-5*p.se,p.hat+5*p.se,length=100)
R> p.samp <- dnorm(pvals,mean=p.hat,sd=p.se)
R> plot(pvals,p.samp,type="l",xlab="",ylab="",
        xlim=p.hat+c(-4,4)*p.se,ylim=c(0,max(p.samp)))
R> abline(h=0,col="gray")
Figure 17-2 gives the result.
Figure 17-2: Visualizing the sampling distribution for the
voting example as per Equation (17.1). The shaded area
represents Pr(0.7 < P^ < 0.75), which is the probability
that the true sample proportion for samples of size
n = 118 lies between 0.7 and 0.75.
Now you can use this distribution to describe the variability in the
sample proportion of voters who already know how they will vote, for
other samples of this size.
For example, the shaded area in Figure 17-2 highlights the probability
that in another sample of the same size, the sample proportion of voters
in the given city who already know how they're going to vote is somewhere
between 0.7 and 0.75. This shaded area can be added with the following code:
  R> pvals.sub <- pvals[pvals>=0.7 & pvals<=0.75]
R> p.samp.sub <- p.samp[pvals>=0.7 & pvals<=0.75]
R> polygon(cbind(c(0.7,pvals.sub,0.75),c(0,p.samp.sub,0)),
           border=NA,col="gray")
Sampling Distributions and Confidence 375
And with knowledge of pnorm, introduced in Section 16.2.2, you can use
the following code to calculate the probability of interest:
  R> pnorm(0.75,mean=p.hat,sd=p.se) - pnorm(0.7,mean=p.hat,sd=p.se)
[1] 0.257238
This sampling distribution suggests that the chance of another sample
proportion, based on the same sample size, lying somewhere between these
two values is about 25.7 percent.
Exercise 17.1
A teacher wants to test all of the 10th-grade students at his school to
gauge their basic mathematical understanding, but the photocopier
breaks after making only six copies of the test. With no other choice,
he chooses six students at random to take the test. Their results,
recorded as a score out of 65, have a sample mean of 41.1. The
standard deviation of the marks of this test is known to be 11.3.
a. Find the standard error associated with the mean test score.
b. Assuming the scores themselves are normally distributed, evaluate the probability that the mean score lies between 45 and 55 if
the teacher took another sample of the same size.
c. A student who gets less than half the questions correct receives a
failing grade (F). Find the probability that the average score is an
F based on another sample of the same size.
A marketing company wants to find out which of two energy drinks
teenagers prefer-drink A or drink B. It surveys 140 teens, and the
results indicate that only 35 percent prefer drink A.
d. Use a quick check to decide whether it is valid to use the normal distribution to represent the sampling distribution of this
proportion.
e. What is the probability that in another sample of the same
size, the proportion of teenagers who prefer drink A is greater
than 0.4?
  f. Find the two values of this sampling distribution that identify the
central 80 percent of values of the proportion of interest.
In Section 16.2.4, the time between cars passing an individual's
location was modeled using an exponential distribution. Say that on
the other side of town, her friend is curious about a similar problem.
Standing outside her house, she records 63 individual times between
cars passing. These sampled times have a mean of x¯ = 37.8 seconds
with a standard deviation of s = 34.51 seconds.
376 Chapter 17
g. The friend inspects a histogram of her raw measurements
and notices that her raw data are heavily right-skewed. Briefly
identify and describe the nature of the sampling distribution
with respect to the sample mean and calculate the appropriate
standard error.
h. Using the standard error from (g) and the appropriate probability distribution, calculate the probability that in another sample
of the same size, the sample mean time between cars passing is as
follows:
  i. More than 40 seconds
ii. Less than half a minute
iii. Between the given sample mean and 40 seconds
17.1.3 Sampling Distributions for Other Statistics
So far you've looked at sampling distributions in cases dealing with a single
sample mean or sample proportion, though it's important to note that many
problems require more complicated measures. Nevertheless, you can apply
the ideas explored in this section to any statistic estimated from a finite-sized
sample. The key, always, is to be able to understand the variability associated
with your point estimates.
In some settings, such as those covered so far, the sampling distribution is parametric, meaning that the functional (mathematical) form of the
probability distribution itself is known and depends only on the provision of
specific parameter values. This is sometimes contingent upon the satisfaction of certain conditions, as you've seen with the application of the normal
distribution covered in this chapter. For other statistics, it may be the case
that you do not know the form of the appropriate sampling distribution-
in these cases, you could use computer simulation to obtain the required
probabilities.
In the remainder of this chapter and over the next few chapters, you'll
continue to explore statistics that are tied to parametric sampling distributions for common tests and models.
NOTE The variability of an estimated quantity is actually only one side of the coin. Just as
important is the issue of statistical bias. Where "natural variability" should be associated with random error, bias is associated with systematic error, in the sense that
a biased statistic does not settle on the corresponding true parameter value as the sample size increases. Bias can be caused by flaws in a study design or collection of data or
can be the result of a poor estimator of the statistic of interest. Bias is an undesirable
trait of any given estimator and/or statistical analysis unless it can be quantified and
removed, which is often difficult if not impossible in practice. I've therefore dealt so far
only with unbiased statistical estimators, many of which are those you may already be
familiar with (for example, the arithmetic mean), and I'll continue to assume unbiasedness moving forward.
Sampling Distributions and Confidence 377
17.2 Confidence Intervals
A confidence interval (CI) is an interval defined by a lower limit l and an upper
limit u, used to describe possible values of a corresponding true population parameter in light of observed sample data. Interpretation of a confidence interval therefore allows you to state a "level of confidence" that a
true parameter of interest falls between this upper and lower limit, often
expressed as a percentage. As such, it is a common and useful tool built
directly from the sampling distribution of the statistic of interest.
The following are the important points to note:
  . The level of confidence is usually expressed as a percentage, such that
you'd construct a 100 × (1 ??? ??) percent confidence interval, where
0 < ?? < 1 is an "amount of tail probability."
. The three most common intervals are defined with either ?? = 0.1 (a
                                                                   90 percent interval), ?? = 0.05 (a 95 percent interval), or ?? = 0.01
(a 99 percent interval).
. Colloquially, you'd state the interpretation of a confidence interval (l,u)
as "I am 100 × (1 ??? ??) percent confident that the true parameter value
lies somewhere between l and u."
Confidence intervals may be constructed in different ways, depending
on the type of statistic and therefore the shape of the corresponding sampling distribution. For symmetrically distributed sample statistics, like those
involving means and proportions that will be used in this chapter, a general
formula is
statistic ± critical value × standard error, (17.2)
where statistic is the sample statistic under scrutiny, critical value is a value
from the standardized version of the sampling distribution that corresponds
to ??, and standard error is the standard deviation of the sampling distribution. The product of the critical value and standard error is referred to
as the error component of the interval; subtraction of the error component
from the value of the statistic provides l, and addition provides u.
With reference to the appropriate sampling distribution, all that a CI
yields are the two values of the distribution that mark off the central 100 ×
(1 ??? ??) percent of the area under the density. (This is the process that was
                                                briefly mentioned in Exercise 17.1 (f).) You then use the CI to make further
interpretations concerning the true (typically unknown) parameter value
that's being estimated by the statistic of interest.
17.2.1 An Interval for a Mean
You know from Section 17.1.1 that the sampling distribution of a single
sample mean depends primarily on whether you know the true standard
deviation of the raw measurements, ??X. Then, provided the sample size for
this sample mean is roughly n ??? 30, the CLT ensures a symmetric sampling
distribution-which will be normal if you know the true value of ??X, or t
based with ?? = n ??? 1 df if you must use the sample standard deviation, s, to
378 Chapter 17
estimate ??X (as is more common in practice). You've seen that the standard
error is defined as the standard deviation divided by the square root of n.
For a small n, you must also assume that the raw observations are normally
distributed, since the CLT will not apply.
To construct an appropriate interval, you must first find the critical
value corresponding to ??. By definition the CI is symmetric, so this translates
to a central probability of (1 ??? ??) around the mean, which is exactly ??/2 in
the lower tail and the same in the upper tail.
Return to the example from Section 17.1.1, dealing with the mean daily
maximum temperatures (degrees Celsius) in January for Dunedin, New
Zealand. Suppose you know the observations are normally distributed but
you don't know the true mean µX (which is set at 22) or the true standard
deviation ??X (which is set at 1.5). Setting it up in the same way as earlier,
assume you've made the following five independent observations:
  R> temp.sample <- rnorm(n=5,mean=22,sd=1.5)
R> temp.sample
[1] 20.46097 21.45658 21.06410 20.49367 24.92843
As you're interested in the sample mean and its sampling distribution,
you must calculate the sample mean x¯, the sample standard deviation s, and
the appropriate standard error of the sample mean, s/
  ???
n.
R> temp.mean <- mean(temp.sample)
R> temp.mean
[1] 21.68075
R> temp.sd <- sd(temp.sample)
R> temp.sd
[1] 1.862456
R> temp.se <- temp.sd/sqrt(5)
R> temp.se
[1] 0.8329155
Now, let's say the aim is to construct a 95 percent confidence interval
for the true, unknown mean µX. This implies ?? = 0.05 (the total amount of
                                                      tail probability) for the relevant sampling distribution. Given the fact that
you know the raw observations are normal and that you're using s (not ??X),
the appropriate distribution is the t-distribution with n ??? 1 = 4 degrees of
freedom. For a central area of 0.95 under this curve, ??/2 = 0.025 must be
in either tail. Knowing that R's q functions operate based on a total lower tail
area, the (positive) critical value is therefore found by supplying a probability of 1 ??? ??/2 = 0.975 to the appropriate function.
R> 1-0.05/2
[1] 0.975
R> critval <- qt(0.975,df=4)
R> critval
[1] 2.776445
Sampling Distributions and Confidence 379
Figure 17-3 shows why the qt function is used in this way (since I used
                                                           similar code throughout Chapter 16, I haven't reproduced the code for Figure 17-3 here).
Figure 17-3: Illustrating the role of the critical value in a confidence
interval for a sample mean, using the Dunedin temperature example.
The sampling distribution is t with 4 df, and the use of qt with respect
to symmetric tail probabilities related to ??/2 = 0.025 yields a central
area of 0.95.
Note that when viewed with respect to the negative version of the
same critical value ("reflected" around the mean and obtained by using
                      qt(0.025,4)), the central, symmetric area under the curve must be 0.95. You
can confirm this using pt.
R> pt(critval,4)-pt(-critval,4)
[1] 0.95
So, all the ingredients are present. You find the 95 percent confidence
interval for the true mean µX via Equation (17.2) with the following lines,
which give l and u, respectively:
  R> temp.mean-critval*temp.se
[1] 19.36821
R> temp.mean+critval*temp.se
[1] 23.99329
380 Chapter 17
The CI given by (19.37,23.99) is therefore interpreted as follows: you
are 95 percent confident that the true mean maximum temperature in
Dunedin in January lies somewhere between 19.37 and 23.99 degrees
Celsius.
With this result, you've combined knowledge of the estimate of the
mean itself with the inherent variability of a sample to define an interval of
values in which you're fairly sure the true mean will lie. As you know, the
true mean in this case is 22, which is indeed included in the calculated CI.
From this, it's easy to alter the intervals to change the confidence levels.
You need to change only the critical value, which, as always, must define ??/2
in each tail. For example, an 80 percent CI (?? = 0.2) and a 99 percent CI
(?? = 0.01) for the same example value given here can be found with these
two lines, respectively:
  R> temp.mean+c(-1,1)*qt(p=0.9,df=4)*temp.se
[1] 20.40372 22.95778
R> temp.mean+c(-1,1)*qt(p=0.995,df=4)*temp.se
[1] 17.84593 25.51557
Note here the use of multiplication by the vector c(-1,1) so that the
lower and upper limits can be obtained at once and the result returned as
a vector of length 2. As usual, the qt function is used with respect to a complete lower-tail area, so p is set at 1 ??? ??/2.
These most recent intervals highlight the natural consequence of moving to a higher confidence level for a given CI. A higher probability in the
central area translates directly to a more extreme critical value, resulting in
a wider interval. This makes sense-in order to be "more confident" about
the true parameter value, you'd need to take into account a larger range of
possible values.
17.2.2 An Interval for a Proportion
Establishing a CI for a sample proportion follows the same rules as for the
mean. With knowledge of the sampling distribution as per Section 17.1.2,
you obtain critical values from the standard normal distribution, and for an
estimate of p^ from a sample of size n, the interval itself is constructed with
the standard error p
p^(1 ??? p^)/n.
Let's return to the example from Section 17.1.2, where 80 of 118 surveyed individuals said that they knew how they were going to vote in the next
US presidential election. Recall you have the following:
  R> p.hat <- 80/118
R> p.hat
[1] 0.6779661
R> p.se <- sqrt(p.hat*(1-p.hat)/118)
R> p.se
[1] 0.04301439
Sampling Distributions and Confidence 381
To construct a 90 percent CI (?? = 0.1), the appropriate critical value
from the standardized sampling distribution of interest is as follows, implying Pr(???1.644854 < Z < 1.644854) = 0.9 for Z ??? N(0,1):
  R> qnorm(0.95)
[1] 1.644854
Now you again follow Equation (17.2):
  R> p.hat+c(-1,1)*qnorm(0.95)*p.se
[1] 0.6072137 0.7487185
You can conclude that you're 90 percent confident that the true proportion of voters who know how they will vote in the next election lies somewhere between 0.61 and 0.75 (rounded to two decimal places).
17.2.3 Other Intervals
The two simple situations presented in Sections 17.2.1 and 17.2.2 serve to
highlight the importance of associating any point estimate (in other words,
                                                            a sample statistic) with the idea of its variability. Confidence intervals can of
course be constructed for other quantities, and over the following sections
(as part of testing hypotheses), I'll expand on the discussion of confidence
intervals to investigate differences between two means and two proportions,
as well as ratios of categorical counts. These more complicated statistics
come with their own standard error formulas, though the corresponding
sampling distributions are still symmetric via the normal and t-curves (if,
                                                                        again, some standard assumptions are met), which means that the now
familiar formulation of Equation (17.2) still applies.
Generally, a confidence interval seeks to mark off a central area of 1 ??? ??
from the sampling distribution of interest, including sampling distributions
that are asymmetric. In those cases, however, it doesn't make much sense
to have a symmetric CI based on a single, standardized critical value as per
Equation (17.2). Similarly, you might not know the functional, parametric
form of the sampling distribution and so may not be willing to make any
distributional assumptions, such as symmetry. In these cases, you can take
an alternative path based on the raw quantiles (or estimated raw quantiles;
                                                see Section 13.2.3) of the supposed asymmetric sampling distribution. Using
specific quantile values to mark off identical ??/2 upper- and lower-tail areas
is a valid method that remains sensitive to the shape of the sampling distribution of interest, while still allowing you to construct a useful interval that
describes potential true parameter values.
17.2.4 Comments on Interpretation of a CI
The typical statement about the interpretation of any CI references a degree
of confidence in where the true parameter value lies, but a more formally
382 Chapter 17
correct interpretation should consider and clarify the probabilistic nature
of the construction. Technically, given a 100(1 ??? ??) percent confidence
level, the more accurate interpretation is as follows: over many samples of
the same size and from the same population where a CI, of the same confidence level, is constructed with respect to the same statistic from each
sample, you would expect the true corresponding parameter value to fall
within the limits of 100(1 ??? ??) percent of those intervals.
This comes from the fact that the theory of a sampling distribution
describes the variability in multiple samples, not just the sample that has
been taken. At first glance it may be difficult to fully appreciate the difference between this and the colloquially used "confidence statement," but it
is important to remain aware of the technically correct definition, particularly given that a CI is typically estimated based on only one sample.
Exercise 17.2
A casual runner records the average time it takes him to sprint
100 meters. He completes the dash 34 times under identical conditions and finds that the mean of these is 14.22 seconds. Assume
that he knows the standard deviation of his runs is ??X = 2.9 seconds.
a. Construct and interpret a 90 percent confidence interval for the
true mean time.
b. Repeat (a), but this time, assume that the standard deviation is
not known and that s = 2.9 is estimated from the sample. How, if
at all, does this change the interval?
  In a particular country, the true proportion of citizens who are left
handed or ambidextrous is unknown. A random sample of 400
people is taken, and each individual is asked to identify with one
of three options: right-handed only, left-handed only, or ambidextrous. The results show that 37 selected left-handed and 11 selected
ambidextrous.
c. Calculate and interpret a 99 percent CI for the true proportion
of left-handed-only citizens.
d. Calculate and interpret a 99 percent CI for the true proportion
of citizens who are either left-handed or ambidextrous.
In Section 17.2.4, the technical interpretation of a CI with respect to
its confidence level was described as the proportion of many similar
intervals (that is, when calculated for samples of the same size from
           the same population) that contain the true value of the parameter of
interest.
Sampling Distributions and Confidence 383
e. Your task is to write an example to demonstrate this behavior
of confidence intervals using simulation. To do so, follow these
instructions:
  - Set up a matrix (see Chapter 3) filled with NAs (Chapter 6)
that has 5,000 rows and 3 columns.
- Use skills from Chapter 10 to write a for loop that, at each of
5,000 iterations, generates a random sample of size 300 from
an exponential distribution with rate parameter ??e = 0.1
(Section 16.2.4).
- Evaluate the sample mean and sample standard deviation of
each sample, and use these quantities with the critical values
from the appropriate sampling distribution to calculate a
95 percent CI for the true mean of the distribution.
- Within the for loop, the matrix should now be filled, row
by row, with your results. The first column will contain the
lower limit, the second will contain the upper limit, and
the third column will be a logical value that is TRUE if the
corresponding interval contains the true mean of 1/??e and
that is FALSE otherwise.
- When the loop is completed, compute the proportion of
TRUEs in the third column of the filled matrix. You should
find that this proportion is close to 0.95; this will vary randomly each time you rerun the loop.
f. Create a plot that draws the first 100 of your estimated confidence intervals as separate horizontal lines drawn from l to u,
one on top of another. One way to do this is to first create an
empty plot with preset x- and y-limits (the latter as c(1,100)) and
then progressively add each line using lines with appropriate
coordinates (this could be done using another for loop). As a
final touch, add to the plot a red vertical line that denotes the
true mean. Confidence intervals that do not include the true
mean will not intersect that vertical line.
The following shows an example of this plot:
  384 Chapter 17
18
HYPOTHES IS TEST ING
In this chapter, you'll build on your experience with confidence intervals and sampling distributions to make more formal
statements about the value of a true, unknown
parameter of interest. For this, you'll learn about frequentist hypothesis testing, where a probability from
a relevant sampling distribution is used as evidence against some claim
about the true value. When a probability is used in this way, it is referred
to as a p-value. In this chapter, I talk about interpreting results for relatively
basic statistics, but you can apply the same concepts to statistics arising from
more complicated methods (such as regression modeling in Chapter 19).
18.1 Components of a Hypothesis Test
To give you an example of hypothesis testing, suppose I told you that 7 percent of a certain population was allergic to peanuts. You then randomly
selected 20 individuals from that population and found that 18 of them were
allergic to peanuts. Assuming your sample was unbiased and truly reflective
of the population, what would you then think about my claim that the true
proportion of allergic individuals is 7 percent?
  Naturally, you would doubt the correctness of my claim. In other words,
there is such a small probability of observing 18 or more successes out of
20 trials for a set success rate of 0.07 that you can state that you have statistical evidence against the claim that the true rate is 0.07. Indeed, when
defining X as the number of allergic individuals out of 20 by assuming
X ??? BIN(20,0.07), evaluating Pr(X ??? 18) gives you the precise p-value,
which is tiny.
R> dbinom(18,size=20,prob=0.07) + dbinom(19,size=20,prob=0.07) +
  dbinom(20,size=20,prob=0.07)
[1] 2.69727e-19
This p-value represents the probability of observing the results in your
sample, X = 18, or a more extreme outcome (X = 19 or X = 20), if the
chance of success was truly 7 percent.
Before looking at specific hypothesis tests and their implementation in
R, this section will introduce terminology that you'll come across often in
the reporting of such tests.
18.1.1 Hypotheses
As the name would suggest, in hypothesis testing, formally stating a claim
and the subsequent hypothesis test is done with a null and an alternative hypothesis. The null hypothesis is interpreted as the baseline or nochange hypothesis and is the claim that is assumed to be true. The alternative hypothesis is the conjecture that you're testing for, against the null
hypothesis.
In general, null and alternative hypotheses are denoted H0 and HA,
respectively, and they are written as follows:
  H0 : . . .
HA : . . .
The null hypothesis is often (but not always) defined as an equality, =,
to a null value. Conversely, the alternative hypothesis (the situation you're
                                                         testing for) is often defined in terms of an inequality to the null value.
. When HA is defined in terms of a less-than statement, with <, it is onesided; this is also called a lower-tailed test.
. When HA is defined in terms of a greater-than statement, with >, it is
one-sided; this is also called an upper-tailed test.
. When HA is merely defined in terms of a different-to statement, with ,,
it is two-sided; this is also called a two-tailed test.
These test variants are entirely situation specific and depend upon the
problem at hand.
386 Chapter 18
18.1.2 Test Statistic
Once the hypotheses are formed, sample data are collected, and statistics
are calculated according to the parameters detailed in the hypotheses. The
test statistic is the statistic that's compared to the appropriate standardized
sampling distribution to yield the p-value.
A test statistic is typically a standardized or rescaled version of the
sample statistic of interest. The distribution and extremity (that is, distance from zero) of the test statistic are the sole drivers of the smallness of
the p-value (which indicates the strength of the evidence against the null
             hypothesis-see Section 18.1.3). Specifically, the test statistic is determined
by both the difference between the original sample statistic and the null
value and the standard error of the sample statistic.
18.1.3 p-value
The p-value is the probability value that's used to quantify the amount of
evidence, if any, against the null hypothesis. More formally, the p-value is
found to be the probability of observing the test statistic, or something more
extreme, assuming the null hypothesis is true.
The exact nature of calculating a p-value is dictated by the type of statistics being tested and the nature of HA. In reference to this, you'll see the
following terms:
  . A lower-tailed test implies the p-value is a left-hand tail probability from
the sampling distribution of interest.
. For an upper-tailed test, the p-value is a right-hand tail probability.
. For a two-sided test, the p-value is the sum of a left-hand tail probability and right-hand tail probability. When the sampling distribution is
symmetric (for example, normal or t, as in all examples coming up in
           Sections 18.2 and 18.3), this is equivalent to two times the area in one
of those tails.
Put simply, the more extreme the test statistic, the smaller the p-value.
The smaller the p-value, the greater the amount of statistical evidence
against the assumed truth of H0.
18.1.4 Significance Level
For every hypothesis test, a significance level, denoted ??, is assumed. This is
used to qualify the result of the test. The significance level defines a cutoff
point, at which you decide whether there is sufficient evidence to view H0 as
incorrect and favor HA instead.
. If the p-value is greater than or equal to ??, then you conclude there
is insufficient evidence against the null hypothesis, and therefore you
retain H0 when compared to HA.
Hypothesis Testing 387
. If the p-value is less than ??, then the result of the test is statistically significant. This implies there is sufficient evidence against the null hypothesis,
and therefore you reject H0 in favor of HA.
Common or conventional values of ?? are ?? = 0.1, ?? = 0.05, and
?? = 0.01.
18.1.5 Criticisms of Hypothesis Testing
The terminology just presented becomes easier to understand once you
look at some examples in the upcoming sections. However, even at this early
stage, it's important to recognize that hypothesis testing is susceptible to
justifiable criticism. The end result of any hypothesis test is to either retain
or reject the null hypothesis, a decision that is solely dependent upon the
rather arbitrary choice of significance level ??; this is most often simply set at
one of the conventionally used values.
Before you begin looking at examples, it is also important to note that a
p-value never provides "proof" of either H0 or HA being truly correct. It can
only ever quantify evidence against the null hypothesis, which one rejects
given a sufficiently small p-value < ??. In other words, rejecting a null hypothesis is not the same as disproving it. Rejecting H0 merely implies that the
sample data suggest HA ought to be preferred, and the p-value merely indicates the strength of this preference.
In recent years, there has been a push against emphasizing these aspects
of statistical inference in some introductory statistics courses owing at least
in part to the overuse, and even misuse, of p-values in some areas of applied
research. A particularly good article by Sterne and Smith (2001) discusses
the role of, and problems surrounding, hypothesis testing from the point of
view of medical research. Another good reference is Reinhart (2015), which
discusses common misinterpretations of p-values in statistics.
That being said, probabilistic inference with respect to sampling distributions is, and will always remain, a cornerstone of frequentist statistical
practice. The best way to improve the use and interpretation of statistical
tests and modeling is with a sound introduction to the relevant ideas and
methods so that, from the outset, you understand statistical significance and
what it can and cannot tell you.
18.2 Testing Means
The validity of hypothesis tests involving sample means is dependent upon
the same assumptions and conditions mentioned in Section 17.1.1. In particular, throughout this section you should assume that the central limit
theorem holds, and if the sample sizes are small (in other words, roughly
                                                  less than 30), the raw data are normally distributed. You'll also focus on
examples where the sample standard deviation s is used to estimate the
true standard deviation, ??X, because this is the most common situation
you'll encounter in practice. Again, mirroring Section 17.1.1, this means
you need to use the t-distribution instead of the normal distribution when
calculating the critical values and p-values.
388 Chapter 18
18.2.1 Single Mean
As you've already met the standard error formula, s/
  ???
n, and the R functionality needed to obtain quantiles and probabilities from the t-distribution (qt
                                                                                                 and pt), the only new concepts to introduce here are related to the definition of the hypotheses themselves and the interpretation of the result.
Calculation: One-Sample t-Test
Let's dive straight into an example-a one-sample t-test. Recall the problem in Section 16.2.2 where a manufacturer of a snack was interested in the
mean net weight of contents in an advertised 80-gram pack. Say that a consumer calls in with a complaint-over time they have bought and precisely
weighed the contents of 44 randomly selected 80-gram packs from different
stores and recorded the weights as follows:
  R> snacks <- c(87.7,80.01,77.28,78.76,81.52,74.2,80.71,79.5,77.87,81.94,80.7,
                 82.32,75.78,80.19,83.91,79.4,77.52,77.62,81.4,74.89,82.95,
                 73.59,77.92,77.18,79.83,81.23,79.28,78.44,79.01,80.47,76.23,
                 78.89,77.14,69.94,78.54,79.7,82.45,77.29,75.52,77.21,75.99,
                 81.94,80.41,77.7)
The customer claims that they've been shortchanged because their data
cannot have arisen from a distribution with mean µ = 80, so the true mean
weight must be less than 80. To investigate this claim, the manufacturer conducts a hypothesis test using a significance level of ?? = 0.05.
First, the hypotheses must be defined, with a null value of 80 grams.
Remember, the alternative hypothesis is "what you're testing for"; in this
case, HA is that µ is smaller than 80. The null hypothesis, interpreted as "no
change," will be defined as µ = 80: that the true mean is in fact 80 grams.
These hypotheses are formalized like this:
  H0 : µ = 80
HA : µ < 80 (18.1)
Second, the mean and standard deviation must be estimated from the
sample.
R> n <- length(snacks)
R> snack.mean <- mean(snacks)
R> snack.mean
[1] 78.91068
R> snack.sd <- sd(snacks)
R> snack.sd
[1] 3.056023
The question your hypotheses seek to answer is this: given the estimated
standard deviation, what's the probability of observing a sample mean (when
                                                                       n = 44) of 78.91 grams or less if the true mean is 80 grams? To answer this,
you need to calculate the relevant test statistic.
Hypothesis Testing 389
Formally, the test statistic T in a hypothesis test for a single mean with
respect to a null value of µ0 is given as
T =
  x¯ ??? µ0
(s/
    ???
  n)
(18.2)
based on a sample of size n, a sample mean of x¯, and a sample standard deviation of s (the denominator is the estimated standard error of the mean).
Assuming the relevant conditions have been met, T follows a t-distribution
with ?? = n ??? 1 degrees of freedom.
In R, the following provides you with the standard error of the sample
mean for the snacks data:
  R> snack.se <- snack.sd/sqrt(n)
R> snack.se
[1] 0.4607128
Then, T can be calculated as follows:
  R> snack.T <- (snack.mean-80)/snack.se
R> snack.T
[1] -2.364419
Finally, the test statistic is used to obtain the p-value. Recall that the
p-value is the probability that you observe T or something more extreme.
The nature of "more extreme" is determined by the alternative hypothesis
HA, which, as a less-than statement, directs you to find a left-hand, lower-tail
probability as the p-value. In other words, the p-value is provided as the area
under the sampling distribution (a t-distribution with 43 df in the current
                                 example) to the left of a vertical line at T. From Section 16.2.3, this is easily
done, as shown here:
  R> pt(snack.T,df=n-1)
[1] 0.01132175
Your result states that if the H0 were true, there would be only a little
more than a 1 percent chance that you'd observe the customer's sample
mean of x¯ = 78.91, or less, as a random phenomenon. Since this p-value is
smaller than the predefined significance level of ?? = 0.05, the manufacturer concludes that there is sufficient evidence to reject the null hypothesis
in favor of the alternative, suggesting the true value of µ is in fact less than
80 grams.
Note that if you find the corresponding 95 percent CI for the single
sample mean, as described in Section 17.2.1 and given by
R> snack.mean+c(-1,1)*qt(0.975,n-1)*snack.se
[1] 77.98157 79.83980
390 Chapter 18
it does not include the null value of 80, mirroring the result of the hypothesis test at the 0.05 level.
R Function: t.test
The result of the one-sample t-test can also be found with the built-in t.test
function.
R> t.test(x=snacks,mu=80,alternative="less")
One Sample t-test
data: snacks
t = -2.3644, df = 43, p-value = 0.01132
alternative hypothesis: true mean is less than 80
95 percent confidence interval:
  -Inf 79.68517
sample estimates:
  mean of x
78.91068
The function takes the raw data vector as x, the null value for the mean
as mu, and the direction of the test (in other words, how to find the p-value
                                      under the appropriate t-curve) as alternative. The alternative argument
has three available options: "less" for HA with <; "greater" for HA with >;
and "two.sided" for HA with ,. The default value of ?? is 0.05. If you want a
different significance level than 0.05, this must be provided to t.test as 1??? ??,
passed to the argument conf.level.
Note that the value of T is reported in the output of t.test, as are the
degrees of freedom and the p-value. You also get a 95 percent "interval,"
but its values of -Inf and 79.68517 do not match the interval calculated
just a moment ago. The manually calculated interval is in fact a two-sided
interval-a bounded interval formed by using an error component that's
equal on both sides.
The CI in the t.test output, on the other hand, takes instruction from
the alternative argument. It provides a one-sided confidence bound. For a lowertailed test, it provides an upper bound on the statistic such that the entire
lower-tail area of the sampling distribution of interest is 0.95, as opposed to
a central area as the traditional two-sided interval does. One-sided bounds
are less frequently used than the fully bounded two-sided interval, which can
be obtained (as the component conf.int) from a relevant call to t.test by
setting alternative="two.sided".
R> t.test(x=snacks,mu=80,alternative="two.sided")$conf.int
[1] 77.98157 79.83980
attr(,"conf.level")
[1] 0.95
Hypothesis Testing 391
This result matches your manually computed version from earlier. Note
also that the corresponding confidence level, 1 ??? ??, is stored alongside this
component as an attribute (refer to Section 6.2.1).
In examining the result for the snack example, with a p-value of around
0.011, remember to be careful when interpreting hypothesis tests. With ?? set
at 0.05 for this particular test, H0 is rejected. But what if the test were carried
out with ?? = 0.01? The p-value is greater than 0.01, so in that case, H0 would
be retained, for no reason other than the arbitrary movement of the value of
??. In these situations, it's helpful to comment on the perceived strength of
the evidence against the null hypothesis. For the current example, you could
reasonably state that there exists some evidence to support HA but that this
evidence is not especially strong.
Exercise 18.1
a. Adult domestic cats of a certain breed are said to have an average
weight of 3.5 kilograms. A feline enthusiast disagrees and collects
a sample of 73 weights of cats of this breed. From her sample,
she calculates a mean of 3.97 kilograms and a standard deviation
of 2.21 kilograms. Perform a hypothesis test to test her claim
that the true mean weight µ is not 3.5 kilograms by setting up the
appropriate hypothesis, carrying out the analysis, and interpreting the p-value (assume the significance level is ?? = 0.05).
b. Suppose it was previously believed that the mean magnitude of
seismic events off the coast of Fiji is 4.3 on the Richter scale. Use
the data in the mag variable of the ready-to-use quakes data set,
providing 1,000 sampled seismic events in that area, to test the
claim that the true mean magnitude is in fact greater than 4.3.
Set up appropriate hypotheses, use t.test (conduct the test at a
                                           significance level of ?? = 0.01), and draw a conclusion.
c. Manually compute a two-sided confidence interval for the true
mean of (b).
18.2.2 Two Means
Often, testing a single sample mean isn't enough to answer the question
you're interested in. In many settings, a researcher wants to directly compare the means of two distinct groups of measurements, which boils down
to a hypothesis test for the true difference between two means; call them µ1
and µ2.
The way in which two groups of data relate to each other affects the specific form of standard error for the difference between two sample means
and therefore the test statistic itself. The actual comparison of the two
392 Chapter 18
means, however, is often of the same nature-the typical null hypothesis is
usually defined as µ1 and µ2 being equal. In other words, the null value of
the difference between the two means is often zero.
Unpaired/Independent Samples: Unpooled Variances
The most general case is where the two sets of measurements are based on
two independent, separate groups (also referred to as unpaired samples).
You compute the sample means and sample standard deviations of both data
sets, define the hypotheses of interest, and then calculate the test statistic.
When you cannot assume the variances of the two populations are
equal, then you perform the unpooled version of the two-sample t-test; this
will be discussed first. If, however, you can safely assume equal variances,
then you can perform a pooled two-sample t-test, which improves the precision of the results. You'll look at the pooled version of the test in a moment.
For an unpooled example, return to the 80-gram snack packet example
from Section 18.2.1. After collecting a sample of 44 packs from the original
manufacturer (label this sample size n1), the disgruntled consumer goes out
and collects n2 = 31 randomly selected 80-gram packs from a rival snack
manufacturer. This second set of measurements is stored as snacks2.
R> snacks2 <- c(80.22,79.73,81.1,78.76,82.03,81.66,80.97,81.32,80.12,78.98,
                79.21,81.48,79.86,81.06,77.96,80.73,80.34,80.01,81.82,79.3,
                79.08,79.47,78.98,80.87,82.24,77.22,80.03,79.2,80.95,79.17,81)
From Section 18.2.1, you already know the mean and standard deviation
of the first sample of size n1 = 44-these are stored as snack.mean (around
                                                                    78.91) and snack.sd (around 3.06), respectively-think of these as x¯1 and s1.
Compute the same quantities, x¯2 and s2, respectively, for the new data.
R> snack2.mean <- mean(snacks2)
R> snack2.mean
[1] 80.1571
R> snack2.sd <- sd(snacks2)
R> snack2.sd
[1] 1.213695
Let the true mean of the original sample be denoted with µ1 and the
true mean of the new sample from the rival company packs be denoted with
µ2. You're now interested in testing whether there is statistical evidence to
support the claim that µ2 is greater than µ1. This suggests the hypotheses of
H0 : µ1 = µ2 and HA : µ1 < µ2, which can be written as follows:
  H0 : µ2 ??? µ1 = 0
HA : µ2 ??? µ1 > 0
That is, the difference between the true mean of the rival company
packs and the original manufacturer's packs, when the original is subtracted
Hypothesis Testing 393
from the rival, is bigger than zero. The "no-change" scenario, the null
hypothesis, is that the two means are the same, so their difference is
truly zero.
Now that you've constructed the hypotheses, let's look at how to actually test them. The difference between two means is the quantity of interest.
For two independent samples arising from populations with true means µ1
and µ2, sample means x¯1 and x¯2, and sample standard deviations s1 and s2,
respectively (and that meet the relevant conditions for the validity of the
              t-distribution), the standardized test statistic T for testing the difference
between µ2 and µ1, in that order, is given as
T =
  x¯2 ??? x¯1 ??? µ0
q
s
2
1
/n1 + s
2
2
/n2
, (18.3)
whose distribution is approximated by a t-distribution with ?? degrees of freedom, where
?? =
  $
  (s
   2
   1
   /n1 + s
   2
   2
   /n2)
2
(s
  2
  1
  /n1)
2/(n1 ??? 1) + (s
              2
              2
              /n2)
2/(n2 ??? 1)
%
(18.4)
In (18.3), µ0 is the null value of interest-typically zero in tests concerned with "difference" statistics. This term would therefore disappear
from the numerator of the test statistic. The denominator of T is the standard error of the difference between two means in this setting.
The ??? · ??? on the right of Equation (18.4) denotes a floor operation-
rounding strictly down to the nearest integer.
NOTE This two-sample t-test, conducted using Equation (18.3), is also called Welch's
t-test. This refers to use of Equation (18.4), called the Welch-Satterthwaite equation. Crucially, it assumes that the two samples have different true variances, which
is why it's called the unpooled variance version of the test.
It's important to be consistent when defining your two sets of parameters and when constructing the hypotheses. In this example, since the test
aims to find evidence for µ2 being greater than µ1, a difference of µ2???µ1 > 0
forms HA (a greater-than, upper-tailed test), and this order of subtraction
is mirrored when calculating T. The same test could be carried out if you
defined the difference the other way around. In that case, your alternative
hypothesis would suggest a lower-tailed test because if you're testing for
µ2 being bigger than µ1, HA would correctly be written as µ1 ??? µ2 < 0.
Again, this would modify the order of subtraction in the numerator of Equation (18.3) accordingly.
The same care must apply to the use of t.test for two-sample comparisons. The two samples must be supplied as the arguments x and y, but the
function interprets x as greater than y when doing an upper-tailed test and
interprets x as less than y when doing a lower-tailed test. Therefore, when
394 Chapter 18
performing the test with alternative="greater" for the snack pack example,
it's snacks2 that must be supplied to x:
  R> t.test(x=snacks2,y=snacks,alternative="greater",conf.level=0.9)
Welch Two Sample t-test
data: snacks2 and snacks
t = 2.4455, df = 60.091, p-value = 0.008706
alternative hypothesis: true difference in means is greater than 0
90 percent confidence interval:
  0.5859714 Inf
sample estimates:
  mean of x mean of y
80.15710 78.91068
With a small p-value of 0.008706, you'd conclude that there is sufficient
evidence to reject H0 in favor of HA (indeed, the p-value is certainly smaller
                                      than the stipulated ?? = 0.1 significance level as implied by conf.level=0.9).
The evidence suggests that the mean net weight of snacks from the rival
manufacturer's 80-gram packs is greater than the mean net weight for the
original manufacturer.
Note that the output from t.test has reported a df value of 60.091,
which is the unfloored result of (18.4). You also receive a one-sided confidence bound (based on the aforementioned confidence level), triggered by
the one-sided nature of this test. Again, the more common two-sided 90 percent interval is also useful; knowing that ?? = ???60.091??? = 60 and using the
statistic and the standard error of interest (numerator and denominator of
                                              Equation (18.3), respectively), you can calculate it.
R> (snack2.mean-snack.mean) +
  c(-1,1)*qt(0.95,df=60)*sqrt(snack.sd^2/44+snack2.sd^2/31)
[1] 0.3949179 2.0979120
Here, you've used the previously stored sample statistics snack.mean,
snack.sd (the mean and standard deviation of the 44 raw measurements
          from the original manufacturer's sample), snack2.mean, and snack2.sd (the
                                                                                same quantities for the 31 observations corresponding to the rival manufacturer). Note that the CI takes the same form as detailed by Equation (17.2)
on page 378 and that to provide the correct 1 ??? ?? central area, the q-function
for the appropriate t-distribution requires 1 ??? ??/2 as its supplied probability
value. You can interpret this as being "90 percent confident that the true difference in mean net weight between the rival and the original manufacturer
(in that order) is somewhere between 0.395 and 2.098 grams." The fact that
zero isn't included in the interval, and that the interval is wholly positive,
supports the conclusion from the hypothesis test.
Hypothesis Testing 395
Unpaired/Independent Samples: Pooled Variance
In the unpooled variance example just passed, there was no assumption that
the variances of the two populations whose means were being compared
were equal. This is an important note to make because it leads to the use of
(18.3) for the test statistic calculation and (18.4) for the associated degrees
of freedom in the corresponding t-distribution. However, if you can assume
equivalence of variances, the precision of the test is improved-you use a
different formula for the standard error of the difference and for calculating
the associated df.
Again, the quantity of interest is the difference between two means, written as µ2 ??? µ1. Assume you have two independent samples of sizes n1 and n2
arising from populations with true means µ1 and µ2, sample means x¯1 and
x¯2, and sample standard deviations s1 and s2, respectively, and assume that
the relevant conditions for the validity of the t-distribution have been met.
Additionally, assume that the true variances of the samples, ??
2
1
and ??
2
2
, are
equal such that ??
2
p = ??
2
1
= ??
2
2
.
NOTE There is a simple rule of thumb to check the validity of the "equal variance" assumption. If the ratio of the larger sample standard deviation to the smaller sample standard deviation is less than 2, you can assume equal variances. For example, if
s1 > s2, then if s1
s2
< 2, you can use the pooled variance test statistic that follows.
The standardized test statistic T for this scenario is given as
T =
  x¯2 ??? x¯1 ??? µ0
q
s
2
p (1/n1 + 1/n2)
, (18.5)
whose distribution is a t-distribution with ?? = n1 + n2 ??? 2 degrees of freedom,
where
s
2
p =
  (n1 ??? 1)s
2
1
+ (n2 ??? 1)s
2
2
n1 + n2 ??? 2
(18.6)
is the pooled estimate of the variance of all the raw measurements. This is substituted in place of s1 and s2 in the denominator of Equation (18.3), resulting
in Equation (18.5).
All other aspects of the two-sample t-test remain as earlier, including the
construction of appropriate hypotheses, the typical null value of µ0, and the
calculation and interpretation of the p-value.
For the comparison of the two means in the snack pack example, you'd
find it difficult to justify using the pooled version of the t-test. Applying
the rule of thumb, the two estimated standard deviations (s1 ??? 3.06 and
                                                          s2 ??? 1.21 for the original and rival manufacturer's samples, respectively)
have a large-to-small ratio that is greater than 2.
R> snack.sd/snack2.sd
[1] 2.51795
396 Chapter 18
Though this is rather informal, if the assumption cannot reasonably be
made, it's best to stick with the unpooled version of the test.
To illustrate this, let's consider a new example. The intelligence quotient (IQ) is a quantity commonly used to measure how clever a person
is. IQ scores are reasonably assumed to be normally distributed, and the
average IQ of the population is said to be 100. Say that you're interested in
assessing whether there is a difference in mean IQ scores between men and
women, suggesting the following hypotheses where you have nmen = 12 and
nwomen = 20:
  H0 : µmen ??? µwomen = 0
HA : µmen ??? µwomen , 0
You randomly sample the following data:
  R> men <- c(102,87,101,96,107,101,91,85,108,67,85,82)
R> women <- c(73,81,111,109,143,95,92,120,93,89,119,79,90,126,62,92,77,106,
              105,111)
As usual, let's calculate the basic statistics required.
R> mean(men)
[1] 92.66667
R> sd(men)
[1] 12.0705
R> mean(women)
[1] 98.65
R> sd(women)
[1] 19.94802
These give the sample averages x¯men and x¯women, as well as their respective sample standard deviations smen and swomen. Enter the following to
quickly check the ratio of the standard deviations:
  R> sd(women)/sd(men)
[1] 1.652626
You can see that the ratio of the larger sample standard deviation to the
smaller is less than 2, so you could assume equal variances in carrying out
the hypothesis test.
The t.test command also enables the pooled two-sample t-test as per
Equations (18.5) and (18.6). To execute it, you provide the optional argument var.equal=TRUE (as opposed to the default var.equal=FALSE, which triggers Welch's t-test).
R> t.test(x=men,y=women,alternative="two.sided",conf.level=0.95,var.equal=TRUE)
Two Sample t-test
Hypothesis Testing 397
data: men and women
t = -0.9376, df = 30, p-value = 0.3559
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -19.016393 7.049727
sample estimates:
  mean of x mean of y
92.66667 98.65000
Note also that HA for this example implies a two-tailed test, hence the
provision of alternative="two.sided".
The resulting p-value of this test, 0.3559, is certainly larger than the
conventional cutoff level of 0.05. Thus, your conclusion here is that there
is no evidence to reject H0-there is insufficient evidence to support a true
difference in the mean IQ scores of men compared to women.
Paired/Dependent Samples
Finally, we'll look comparing two means in paired data. This setting is distinctly different from that of both unpaired t-tests because it concerns the
way the data have been collected. The issue concerns dependence between
pairs of observations across the two groups of interest-previously, the measurements in each group have been defined as independent. This notion
has important consequences for how the test can be carried out.
Paired data occur if the measurements forming the two sets of observations are recorded on the same individual or if they are related in some
other important or obvious way. A classic example of this is "before" and
"after" observations, such as two measurements made on each person before
and after some kind of intervention treatment. These situations still focus on
the difference between the mean outcomes in each group, but rather than
working with the two data sets separately, a paired t-test works with a single
mean-the true mean of the individual paired differences µd.
As an example, consider a company interested in the efficacy of a drug
designed to reduce resting heart rates in beats per minute (bpm). The resting heart rates of 16 individuals are measured. The individuals are then
administered a course of the treatment, and their resting heart rates are
again measured. The data are provided in the two vectors rate.before and
rate.after as follows:
  R> rate.before <- c(52,66,89,87,89,72,66,65,49,62,70,52,75,63,65,61)
R> rate.after <- c(51,66,71,73,70,68,60,51,40,57,65,53,64,56,60,59)
It quickly becomes clear why any test comparing these two groups must
take dependence into account. Heart rate is affected by an individual's age,
build, and level of physical fitness. An unfit individual older than 60 is likely
to have a higher baseline resting heart rate than a fit 20-year-old, and if both
are given the same drug to lower their heart rate, their final heart rates are
398 Chapter 18
still likely to reflect the baselines. Any true effect of the drug therefore has
the potential to be hidden if you approached the analysis using either of the
unpaired t-tests.
To overcome this problem, the paired two-sample t-test considers the
difference between each pair of values. Labeling one set of n measurements
as x1, . . ., xn and the other set of n observations as y1, . . ., yn, the difference,
d, is defined as di = yi ??? xi
; i = 1, . . ., n. In R, you can easily compute the
pairwise differences:
  R> rate.d <- rate.after-rate.before
R> rate.d
[1] -1 0 -18 -14 -19 -4 -6 -14 -9 -5 -5 1 -11 -7 -5 -2
The following code calculates the sample mean d¯ and standard deviation sd of these differences:
  R> rate.dbar <- mean(rate.d)
R> rate.dbar
[1] -7.4375
R> rate.sd <- sd(rate.d)
R> rate.sd
[1] 6.196437
You want to see how much the heart rate is reduced by, so the test at
hand will be concerned with the following hypotheses:
  H0 : µd = 0
HA : µd < 0
Given the order or subtraction used to obtain the differences, detection
of a successful reduction in heart rate will be represented by an "after" mean
that is smaller than the "before" mean.
Expressing all this mathematically, the value of interest is the true mean
difference, µd, between two means of dependent pairs of measurements.
There are two sets of n measurements, x1, . . ., xn and y1, . . ., yn, with pairwise differences d1, . . ., dn. The relevant conditions for the validity of the
t-distribution must be met; in this case, if the number of pairs n is less than
30, then you must be able to assume the raw data are normally distributed.
The test statistic T is given as
T =
  d¯??? µ0
sd/
  ???
n
, (18.7)
where d¯ is the mean of the pairwise differences, sd is the sample standard
deviation of the pairwise differences, and µ0 is the null value (usually zero).
The statistic T follows a t-distribution with n ??? 1 df.
The form of Equation (18.7) is actually the same as the form of the
test statistic in (18.2), once the sample statistics for the individual paired
Hypothesis Testing 399
differences have been calculated. Furthermore, it's important to note that
n represents the total number of pairs, not the total number of individual
observations.
For the current example hypotheses, you can find the test statistic and
p-value with rate.dbar and rate.sd.
R> rate.T <- rate.dbar/(rate.sd/sqrt(16))
R> rate.T
[1] -4.801146
R> pt(rate.T,df=15)
[1] 0.000116681
These results suggest evidence to reject H0. In t.test, the optional logical argument paired must be set to TRUE.
R> t.test(x=rate.after,y=rate.before,alternative="less",conf.level=0.95,
          paired=TRUE)
Paired t-test
data: rate.after and rate.before
t = -4.8011, df = 15, p-value = 0.0001167
alternative hypothesis: true difference in means is less than 0
95 percent confidence interval:
  -Inf -4.721833
sample estimates:
  mean of the differences
-7.4375
Note that the order you supply your data vectors to the x and y arguments follows the same rules as for the unpaired tests, given the desired
value of alternative. The same p-value as was calculated manually is confirmed through the use of t.test, and since this is less than an assumed conventional significance level of ?? = 0.05, a valid conclusion would be to state
that there is statistical evidence that the medication does reduce the mean
resting heart rate. You could go on to say you're 95 percent confident that
the true mean difference in heart rate after taking the course of medication
lies somewhere between
R> rate.dbar-qt(0.975,df=15)*(rate.sd/sqrt(16))
[1] -10.73935
and
R> rate.dbar+qt(0.975,df=15)*(rate.sd/sqrt(16))
[1] -4.135652
400 Chapter 18
NOTE On some occasions, such as when your data strongly indicate non-normality, you may
not be comfortable assuming the validity of the CLT (refer back to Section 17.1.1). An
alternative approach to the tests discussed here is to employ a nonparametric technique that relaxes these distributional requirements. In the two-sample case, you could
employ the Mann-Whitney U test (also known as the Wilcoxon rank-sum test).
This is a hypothesis test that compares two medians, as opposed to two means. You
can use the R function wilcox.test to access this methodology; its help page provides
useful commentary and references on the particulars of the technique.
Exercise 18.2
In the package MASS you'll find the data set anorexia, which contains
data on pre- and post-treatment weights (in pounds) of 72 young
women suffering from the disease, obtained from Hand et al. (1994).
One group of women is the control group (in other words, no intervention), and the other two groups are the cognitive behavioral
program and family support intervention program groups. Load the
library and ensure you can access the data frame and understand its
contents. Let µd denote the mean difference in weight, computed as
(post-weight ??? pre-weight).
a. Regardless of which treatment group the participants fall
into, conduct and conclude an appropriate hypothesis test
with ?? = 0.05 for the entire set of weights for the following
hypotheses:
  H0 : µd = 0
HA : µd > 0
b. Next, conduct three separate hypothesis tests using the same
defined hypotheses, based on which treatment group the participants fall into. What do you notice?
  Another ready-to-use data set in R is PlantGrowth (Dobson, 1983),
which records a continuous measure of the yields of a certain plant,
looking at the potential effect of two supplements administered
during growth to increase the yield when compared to a control
group with no supplement.
c. Set up hypotheses to test whether the mean yield for the control
group is less than the mean yield from a plant given either of the
treatments. Determine whether this test should proceed using a
pooled estimate of the variance or whether Welch's t-test would
be more appropriate.
d. Conduct the test and make a conclusion (assuming normality of
                                           the raw observations).
Hypothesis Testing 401
As discussed, there is a rule of thumb for deciding whether to use a
pooled estimate of the variance in an unpaired t-test.
e. Your task is to write a wrapper function that calls t.test after
deciding whether it should be executed with var.equal=FALSE
according to the rule of thumb. Use the following guidelines:
  - Your function should take four defined arguments: x and y
with no defaults, to be treated in the same way as the same
arguments in t.test; and var.equal and paired, with defaults
that are the same as the defaults of t.test.
- An ellipsis (Section 9.2.5) should be included to represent
any additional arguments to be passed to t.test.
- Upon execution, the function should determine whether
paired=FALSE.
* If paired is TRUE, then there is no need to proceed with
the check of a pooled variance.
* If paired is FALSE, then the function should determine
the value for var.equal automatically by using the rule of
thumb.
- If the value of var.equal was set automatically, you can assume
it will override any value of this argument initially supplied
by the user.
- Then, call t.test appropriately.
f. Try your new function on all three examples in the text of Section 18.2.2, ensuring you reach identical results.
18.3 Testing Proportions
A focus on means is especially common in statistical modeling and hypothesis testing, and therefore you must also consider sample proportions, interpreted as the mean of a series of n binary trials, in which the results are
success (1) or failure (0). This section focuses on the parametric tests of
proportions, which assume normality of the target sampling distributions
(otherwise referred to as Z-tests).
The general rules regarding the setup and interpretation of hypothesis
tests for sample proportions remain the same as for sample means. In this
introduction to Z-tests, you can consider these as tests regarding the true
value of a single proportion or the difference between two proportions.
18.3.1 Single Proportion
Section 17.1.2 introduced the sampling distribution of a single sample
proportion to be normally distributed, with a mean centered on the
true proportion ?? and with a standard error of ???
??(1 ??? ??)/n. Provided
the trials are independent and that n isn't "too small" and ?? isn't "too
close" to 0 or 1, those formulas are applicable here.
402 Chapter 18
NOTE A rule of thumb to check the latter condition on n and ?? simply involves checking that
np^ and n(1 ??? p^) are both greater than 5, where p^ is the sample estimate of ??.
It's worth noting that the standard error in the case of hypothesis
tests involving proportions is itself dependent upon ??. This is important-
remember that any hypothesis test assumes satisfaction of H0 in the relevant
calculations. In dealing with proportions, that means when computing the
test statistic, the standard error must make use of the null value ??0 rather
than the estimated sample proportion p^.
I'll clarify this in an example. Suppose an individual fond of a particular
fast-food chain notices that he tends to have an upset stomach within a certain amount of time after having his usual lunch. He comes across the website of a blogger who believes that the chance of getting an upset stomach
shortly after eating that particular food is 20 percent. The individual is curious to determine whether his true rate of stomach upset ?? is any different
from the blogger's quoted value and, over time, visits these fast-food outlets
for lunch on n = 29 separate occasions, recording the success (TRUE) or failure (FALSE) of experiencing an upset stomach. This suggests the following
pair of hypotheses:
  H0 : ?? = 0.2
HA : ?? , 0.2
These may be tested according to the general rules discussed in the following sections.
Calculation: One-Sample Z-Test
In testing for the true value of some proportion of success, ??, let p^ be the
sample proportion over n trials, and let the null value be denoted with ??0.
You find the test statistic with the following:
  Z =
  p^ ??? ??0
q
??0 (1?????0)
n
(18.8)
Provided the aforementioned conditions on the size of n and the value
of ?? can be assumed, Z ??? N(0,1).
The denominator of Equation (18.8), the standard error of the proportion, is calculated with respect to the null value ??0, not p^. As mentioned just
a moment ago, this is to satisfy the assumption of "truth" of H0 as the test is
carried out, so it allows interpretation of the resulting p-value as usual. The
standard normal distribution is used to find the p-value with respect to Z;
the direction underneath this curve is governed by the nature of HA just as
before.
Getting back to the fast-food example, suppose these are the observed
data, where 1 is recorded for an upset stomach and is 0 otherwise.
sick <- c(0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1)
Hypothesis Testing 403
The number of successes and probability of success in this sample are as
follows:
  R> sum(sick)
[1] 8
R> p.hat <- mean(sick)
R> p.hat
[1] 0.2758621
A quick check indicates that as per the rule of thumb, the test is reasonable to carry out:
  R> 29*0.2
[1] 5.8
R> 29*0.8
[1] 23.2
Following Equation (18.8), the test statistic Z for this example is as
follows:
  R> Z <- (p.hat-0.2)/sqrt(0.2*0.8/29)
R> Z
[1] 1.021324
The alternative hypothesis is two-sided, so you compute the corresponding p-value as a two-tailed area under the standard normal curve. With a
positive test statistic, this can be evaluated by doubling the upper-tailed area
from Z.
R> 2*(1-pnorm(Z))
[1] 0.3071008
Assume a conventional ?? level of 0.05. The high p-value given as 0.307
suggests the results in the sample of size 29 are not unusual enough, under
the assumption that the null hypothesis is true, to reject H0. There is insufficient evidence to suggest that the proportion of instances of an upset
stomach that this individual experiences is any different from 0.2 as noted
by the blogger.
You can support this conclusion with a confidence interval. At the level
of 95 percent, you calculate the CI:
  R> p.hat+c(-1,1)*qnorm(0.975)*sqrt(p.hat*(1-p.hat)/29)
[1] 0.1131927 0.4385314
This interval easily includes the null value of 0.2.
404 Chapter 18
R Function: prop.test
Once more, R rescues you from tedious step-by-step calculation. The
ready-to-use prop.test function allows you to perform, among other things,
a single sample proportion test. The function actually performs the test
in a slightly different way, using the chi-squared distribution (which will be
                                                                 explored more in Section 18.4). However, the test is equivalent, and the
resulting p-value from prop.test is identical to the one reached using the Zbased test.
To the prop.test function, as used for a single sample test of a proportion, you provide the number of successes observed as x, the total number
of trials as n, and the null value as p. The two further arguments, alternative
(defining the nature of HA) and conf.level (defining 1 ??? ??), are identical
to the same arguments in t.test and have defaults of "two.sided" and 0.95,
respectively. Lastly, it is recommended to explicitly set the optional argument correct=FALSE if your data satisfy the np^ and n(1 ??? p^) rule of thumb.
For the current example, you perform the test with this code:
  R> prop.test(x=sum(sick),n=length(sick),p=0.2,correct=FALSE)
1-sample proportions test without continuity correction
data: sum(sick) out of length(sick), null probability 0.2
X-squared = 1.0431, df = 1, p-value = 0.3071
alternative hypothesis: true p is not equal to 0.2
95 percent confidence interval:
  0.1469876 0.4571713
sample estimates:
  p
0.2758621
The p-value is the same as you got earlier. Note, however, that the
reported CI is not quite the same (the normal-based interval, dependent
                                   upon the CLT). The CI produced by prop.test is referred to as the Wilson
score interval, which takes into account the direct association that a "probability of success" has with the binomial distribution. For simplicity, you'll
continue to work with normal-based intervals when performing hypothesis
tests involving proportions here.
Note also that, just like t.test, any one-sided test performed with
prop.test will provide only a single-limit confidence bound; you'll see this
in the following example.
18.3.2 Two Proportions
With a basic extension to the previous procedure, by way of a modification
to the standard error, you can compare two estimated proportions from
independent populations. As with the difference between two means, you're
often testing whether the two proportions are the same and thus have a difference of zero. Therefore, the typical null value is zero.
Hypothesis Testing 405
For an example, consider a group of students taking a statistics exam. In
this group are n1 = 233 students majoring in psychology, of whom x1 = 180
pass, and n2 = 197 students majoring in geography, of whom 175 pass. Suppose it is claimed that the geography students have a higher pass rate in
statistics than the psychology students.
Representing the true pass rates for psychology students as ??1 and geography students as ??2, this claim can be statistically tested using a pair of
hypotheses defined as follows:
  H0 : ??2 ??? ??1 = 0
HA : ??2 ??? ??1 > 0
Just as with a comparison between two means, it's important to keep
the order of differencing consistent throughout the test calculations. This
example shows an upper-tailed test.
Calculation: Two-Sample Z-Test
In testing for the true difference between two proportions mathematically,
??1 and ??2, let p^1 = x1/n1 be the sample proportion for x1 successes in n1
trials corresponding to ??1, and the same quantities as p^2 = x2/n2 for ??2.
With a null value of the difference denoted ??0, the test statistic is given by
the following:
  Z =
  p^2 ??? p^1 ??? ??0
q
p
???
(1 ??? p
 ???
)

1
n1
+
  1
n2

(18.9)
Provided you can assume to apply the aforementioned conditions for a
proportion with respect to n1, n2 and ??1, ??2, you can treat Z ??? N(0,1).
There is a new quantity, p
???
, present in the denominator of (18.9). This
is a pooled proportion, given as follows:
  p
??? =
  x1 + x2
n1 + n2
(18.10)
As noted, in this kind of test it is common for the null value, the true
difference in proportions, to be set to zero (in other words, ??0 = 0).
The denominator of Equation (18.9) is itself the standard error of
the difference between two proportions as used in a hypothesis test. The
need to use p
???
lies once more in the fact that H0 is assumed to be true.
Using p^1 and p^2 separately in the denominator of (18.9), in the form of
p
p^1(1 ??? p^1)/n1 + p^2(1 ??? p^2)/n2 (the standard error of the difference between
                                   two proportions outside the confines of a hypothesis test), would violate the
assumed "truth" of H0.
So, returning to the statistics exam taken by the psychology and geography students, you can evaluate the required quantities as such:
  R> x1 <- 180
R> n1 <- 233
406 Chapter 18
R> p.hat1 <- x1/n1
R> p.hat1
[1] 0.7725322
R> x2 <- 175
R> n2 <- 197
R> p.hat2 <- x2/n2
R> p.hat2
[1] 0.8883249
The results indicate sample pass rates of around 77.2 percent for the
psychology students and 88.8 percent for the geography students; this is
a difference of roughly 11.6 percent. From examining the values of p^1, n1
and p^2, n2, you can see that the rule of thumb is satisfied for this test; again,
assume a standard significance level of ?? = 0.05.
The pooled proportion p
???
, following (18.10), is as follows:
  R> p.star <- (x1+x2)/(n1+n2)
R> p.star
[1] 0.8255814
With that you calculate the test statistic Z as per Equation (18.9) with
the following:
  R> Z <- (p.hat2-p.hat1)/sqrt(p.star*(1-p.star)*(1/n1+1/n2))
R> Z
[1] 3.152693
In light of the hypotheses, you find the corresponding p-value as a righthand, upper-tail area from Z underneath the standard normal curve as
follows:
  R> 1-pnorm(Z)
[1] 0.0008088606
You observe a p-value that's substantially smaller than ??, so the formal
decision is of course to reject the null hypothesis in favor of the alternative.
The sample data provide sufficient evidence against H0 such that you can
conclude that evidence exists to support the pass rate for geography students
being higher than the pass rate for psychology students.
R Function: prop.test
Once more, R allows you to perform the test with one line of code using
prop.test. For comparisons of two proportions, you pass the number of successes in each group as a vector of length 2 to x and the respective sample
sizes as another vector of length 2 to n. Note that the order of the entries
must reflect the order of alternative if this is one-sided (in other words,
                                                            Hypothesis Testing 407
                                                            here, the proportion that is to be tested as "greater" corresponds to the
                                                            first elements of x and n). Once more, correct is set to FALSE.
R> prop.test(x=c(x2,x1),n=c(n2,n1),alternative="greater",correct=FALSE)
2-sample test for equality of proportions without continuity correction
data: c(x2, x1) out of c(n2, n1)
X-squared = 9.9395, df = 1, p-value = 0.0008089
alternative hypothesis: greater
95 percent confidence interval:
  0.05745804 1.00000000
sample estimates:
  prop 1 prop 2
0.8883249 0.7725322
The p-value is identical to the one generated by the previous series of
calculations, suggesting a rejection of H0. Since prop.test was called as a
one-sided test, the confidence interval returned provides a single bound.
To provide a two-sided CI for the true difference, it makes sense, considering the outcome of the test, to construct this using the separate p^1 and p^2
instead of using the denominator of (18.9) specifically (which assumes truth
                                                         of H0). The "separate-estimate" version of the standard error of the difference between two proportions was given earlier (in the text beneath Equation (18.10)), and a 95 percent CI is therefore calculated with the following:
  R> (p.hat2-p.hat1) +
  c(-1,1)*qnorm(0.975)*sqrt(p.hat1*(1-p.hat1)/n1+p.hat2*(1-p.hat2)/n2)
[1] 0.04628267 0.18530270
With that, you're 95 percent confident that the true difference between
the proportion of geography students passing the exam and the proportion
of psychology students passing the exam lies somewhere between 0.046 and
0.185. Naturally, the interval also reflects the result of the hypothesis test-it
doesn't include the null value of zero and is wholly positive.
Exercise 18.3
An advertisement for a skin cream claims nine out of ten women who
use it would recommend it to a friend. A skeptical salesperson in a
department store believes the true proportion of women users who'd
recommend it, ??, is much smaller than 0.9. She follows up with 89
random customers who had purchased the skin cream and asks if
they would recommend it to others, to which 71 answer yes.
408 Chapter 18
a. Set up an appropriate pair of hypotheses for this test and
determine whether it will be valid to carry out using the normal
distribution.
b. Compute the test statistic and the p-value and state your conclusion for the test using a significance level of ?? = 0.1.
c. Using your estimated sample proportion, construct a two-sided
90 percent confidence interval for the true proportion of women
who would recommend the skin cream.
The political leaders of a particular country are curious as to the
proportion of citizens in two of its states that support the decriminalization of marijuana. A small pilot survey taken by officials reveals
that 97 out of 445 randomly sampled voting-age citizens residing in
state 1 support the decriminalization and that 90 out of 419 votingage citizens residing in state 2 support the same notion.
d. Letting ??1 denote the true proportion of citizens in support of
decriminalization in state 1, and ??2 the same measure in state 2,
conduct and conclude a hypothesis test under a significance level
of ?? = 0.05 with reference to the following hypotheses:
  H0 : ??2 ??? ??1 = 0
HA : ??2 ??? ??1 , 0
e. Compute and interpret a corresponding CI.
Though there is standard, ready-to-use R functionality for the t-test,
at the time of this writing, there is no similar function for the Z-test
(in other words, the normal-based test of proportions described
  here) except in contributed packages.
f. Your task is to write a relatively simple R function, Z.test, that
can perform a one- or two-sample Z-test, using the following
guidelines:
  - The function should take the following arguments: p1 and
n1 (no default) to pose as the estimated proportion and
sample size; p2 and n2 (both defaulting to NULL) that contain
the second sample proportion and sample size in the event
of a two-sample test; p0 (no default) as the null value; and
alternative (default "two.sided") and conf.level (default 0.95),
to be used in the same way as in t.test.
- When conducting a two-sample test, it should be p1
that is tested as being smaller or larger than p2 when
alternative="less" or alternative="greater", the same as in
the use of x and y in t.test.
- The function should perform a one-sample Z-test using p1,
n1, and p0 if either p2 or n2 (or both) is NULL.
Hypothesis Testing 409
- The function should contain a check for the rule of thumb
to ensure the validity of the normal distribution in both
one- and two-sample settings. If this is violated, the function should still complete but should issue an appropriate
warning message (see Section 12.1.1).
- All that need be returned is a list containing the members
Z (test statistic), P (appropriate p-value-this can be determined by alternative; for a two-sided test, determining
                       whether Z is positive or not can help), and CI (two-sided
                                                                       CI with respect to conf.level).
g. Replicate the two examples in the text of Sections 18.3.1 and
18.3.2 using Z.test; ensure you reach identical results.
h. Call Z.test(p1=0.11,n1=10,p0=0.1) to try your warning message in
the one-sample setting.
18.4 Testing Categorical Variables
The normal-based Z-test is particular to data that are binary in nature. To
statistically test claims regarding more general categorical variables, with
more than two distinct levels, you use the ubiquitous chi-squared test. Pronounced kai, "chi" refers to the Greek symbol ?? and is sometimes noted in
shorthand as the ??
2
test.
There are two common variants of the chi-squared test. The first-a chisquared test of distribution, also called a goodness of fit (GOF) test-is used
when assessing the frequencies in the levels of a single categorical variable.
The second-a chi-squared test of independence-is employed when you're
investigating the relationship between frequencies in the levels of two such
variables.
18.4.1 Single Categorical Variable
Like the Z-test, the one-dimensional chi-squared test is also concerned with
comparing proportions but in a setting where there are more than two proportions. A chi-squared test is used when you have k levels (or categories) of
a categorical variable and want to hypothesize about their relative frequencies to find out what proportion of n observations fall into each defined category. In the following examples, it must be assumed that the categories are
mutually exclusive (in other words, an observation cannot take more than one
                    of the possible categories) and exhaustive (in other words, the k categories
                                                                cover all possible outcomes).
I'll illustrate how hypotheses are constructed and introduce the relevant
ideas and methods with the following example. Suppose a researcher in sociology is interested in the dispersion of rates of facial hair in men of his local
city and whether they are uniformly represented in the male population. He
defines a categorical variable with three levels: clean shaven (1), beard only
410 Chapter 18
or moustache only (2), and beard and moustache (3). He collects data on 53
randomly selected men and finds the following outcomes:
  R> hairy <- c(2,3,2,3,2,1,3,3,2,2,3,2,2,2,3,3,3,2,3,2,2,2,1,3,2,2,2,1,2,2,3,
                2,2,2,2,1,2,1,1,1,2,2,2,3,1,2,1,2,1,2,1,3,3)
Now, the research question asks whether the proportions in each category are equally represented. Let ??1, ??2, and ??3 represent the true proportion of men in the city who fall into groups 1, 2, and 3, respectively. You
therefore seek to test these hypotheses:
  H0 : ??1 = ??2 = ??3 =
  1
3
HA : H0 is incorrect
For this test, use a standard significance level of 0.05.
The appearance of the alternative hypothesis is a little different from
what you've seen so far but is an accurate reflection of the interpretation of
a chi-squared goodness of fit test. In these types of problems, H0 is always
that the proportions in each group are equal to the stated values, and HA
is that the data, as a whole, do not match the proportions defined in the
null. The test is conducted assuming the null hypothesis is true, and evidence against the no-change, baseline setting will be represented as a small
p-value.
Calculation: Chi-Squared Test of Distribution
The quantities of interest are the proportion of n observations in each of k
categories, ??1, . . ., ??k , for a single mutually exclusive and exhaustive categorical variable. The null hypothesis defines hypothesized null values for each
proportion; label these respectively as ??0(1)
, . . ., ??0(k)
. The test statistic ??
2
is
given as
??
2 =
  X
k
i=1
(Oi ??? Ei)
2
Ei
, (18.11)
where Oi
is the observed count and Ei
is the expected count in the ith category; i = 1, . . ., k. The Oi are obtained directly from the raw data, and the
expected counts, Ei = n??0(i)
, are merely the product of the overall sample
size n with the respective null proportion for each category. The result of
??
2
follows a chi-squared distribution (explained further momentarily) with
?? = k ??? 1 degrees of freedom. You usually consider the test to be valid based
on an informal rule of thumb stating that at least 80 percent of the expected
counts Ei should be at least 5.
In this type of chi-squared test, it is important to note the following:
  . The term goodness of fit refers to the proximity of the observed data to
the distribution hypothesized in H0.
Hypothesis Testing 411
. Positive extremity of the result of (18.11) provides evidence against H0.
As such, the corresponding p-value is always computed as an upper-tail area.
. As in the current example, a test for uniformity simplifies the null
hypothesis slightly by having equivalent null proportions ??0 = ??0(1) =
  . . . = ??0(k)
.
. A rejected H0 doesn't tell you about the true values of ??i
. It merely suggests that they do not follow H0 specifically.
The chi-squared distribution relies on specification of a degree of freedom, much like the t-distribution. Unlike a t curve, however, a chi-squared
curve is unidirectional in nature, being defined for non-negative values and
with a positive (right-hand) horizontal asymptote (tail going to zero).
It's this unidirectional distribution that leads to p-values being defined
as upper-tail areas only; decisions like one- or two-tailed areas have no relevance in these types of chi-squared tests. To get an idea of what the density functions actually look like, Figure 18-1 shows three particular curves
defined with ?? = 1, ?? = 5, and ?? = 10 degrees of freedom.
Figure 18-1: Three instances of the chi-squared density
function using differing degrees of freedom values. Note
the positive domain of the function and the "flattening"
and "right-extending" behavior as ?? is increased.
This image was produced using the relevant d-function, dchisq, with ??
passed to the argument df.
R> x <- seq(0,20,length=100)
R> plot(x,dchisq(x,df=1),type="l",xlim=c(0,15),ylim=c(0,0.5),ylab="density")
R> lines(x,dchisq(x,df=5),lty=2)
R> lines(x,dchisq(x,df=10),lty=3)
R> abline(h=0,col="gray")
412 Chapter 18
R> abline(v=0,col="gray")
R> legend("topright",legend=c("df=1","df=5","df=10"),lty=1:3)
The current facial hair example is a test for the uniformity of the distribution of frequencies in the three categories. You can obtain the observed
counts and corresponding proportions with table.
R> n <- length(hairy)
R> n
[1] 53
R> hairy.tab <- table(hairy)
R> hairy.tab
hairy
1 2 3
11 28 14
R> hairy.tab/n
hairy
1 2 3
0.2075472 0.5283019 0.2641509
For computation of the test statistic ??
2
, you have the observed counts
Oi
in hairy.tab. The expected count Ei
is a straightforward arithmetic calculation of the total number of observations multiplied by the null proportion 1/3 (the result stored as expected), giving you the same value for each
category.
These, as well as the contribution of each category to the test statistic,
are nicely presented in a matrix constructed with cbind (Section 3.1.2).
R> expected <- 1/3*n
R> expected
[1] 17.66667
R> hairy.matrix <- cbind(1:3,hairy.tab,expected,
                         (hairy.tab-expected)^2/expected)
R> dimnames(hairy.matrix) <- list(c("clean","beard OR mous.",
                                    "beard AND mous."),
                                  c("i","Oi","Ei","(Oi-Ei)^2/Ei"))
R> hairy.matrix
i Oi Ei (Oi-Ei)^2/Ei
clean 1 11 17.66667 2.5157233
beard OR mous. 2 28 17.66667 6.0440252
beard AND mous. 3 14 17.66667 0.7610063
Note that all the expected counts are comfortably greater than 5, which
satisfies the informal rule of thumb mentioned earlier. In terms of R coding, note also that the single number expected is implicitly recycled to match
the length of the other vectors supplied to cbind and that you've used the
dimnames attribute (refer to Section 6.2.1) to annotate the rows and columns.
Hypothesis Testing 413
The test statistic, as per (18.11), is given as the sum of the (Oi ??? Ei)
2
/Ei
contributions in the fourth column of hairy.matrix.
R> X2 <- sum(hairy.matrix[,4])
R> X2
[1] 9.320755
The corresponding p-value is the appropriate upper-tail area from the
chi-squared distribution with ?? = 3 ??? 1 = 2 degrees of freedom.
R> 1-pchisq(X2,df=2)
[1] 0.009462891
This small p-value provides evidence to suggest that the true frequencies
in the defined categories of male facial hair are not uniformly distributed
in a 1/3,1/3,1/3 fashion. Remember that the test result doesn't give you the
true proportions but only suggests that they do not follow those in H0.
R Function: chisq.test
Like t.test and prop.test, R provides a quick-use function for performing a chi-squared GOF test. The chisq.test function takes the vector of
observed frequencies as its first argument x. For the facial hair example,
this simple line therefore provides the same results as found previously:
  R> chisq.test(x=hairy.tab)
Chi-squared test for given probabilities
data: hairy.tab
X-squared = 9.3208, df = 2, p-value = 0.009463
By default, the function performs a test for uniformity, taking the number of categories as the length of the vector supplied to x. However, suppose that the researcher collecting the facial hair data realizes that he was
doing so in November, a month during which many men grow mustaches in
support of "Mo-vember" to raise awareness of men's health. This changes
thoughts on the true rates in terms of his clean-shaven (1), beard-only or
moustache-only (2), and beard and moustache (3) categories. He now wants
to test the following:
  H0 : ??0(1) = 0.25; ??0(2) = 0.5; ??0(3) = 0.25
HA : H0 is incorrect.
If a GOF test of uniformity is not desired, when the "true" rates across
the categories are not all the same, the chisq.test function requires you
to supply the null proportions as a vector of the same length as x to the p
argument. Naturally, each entry in p must correspond to the categories tabulated in x.
414 Chapter 18
R> chisq.test(x=hairy.tab,p=c(0.25,0.5,0.25))
Chi-squared test for given probabilities
data: hairy.tab
X-squared = 0.5094, df = 2, p-value = 0.7751
With a very high p-value, there is no evidence to reject H0 in this scenario. In other words, there is no evidence to suggest that the proportions
hypothesized in H0 are incorrect.
18.4.2 Two Categorical Variables
The chi-squared test can also apply to the situation in which you have two
mutually exclusive and exhaustive categorical variables at hand-call them
variable A and variable B. It is used to detect whether there might be some
influential relationship (in other words, dependence) between A and B by
looking at the way in which the distribution of frequencies change together
with respect to their categories. If there is no relationship, the distribution
of frequencies in variable A will have nothing to do with the distribution of
frequencies in variable B. As such, this particular variant of the chi-squared
test is called a test of independence and is always performed with the following
hypotheses:
  H0 : Variables A and B are independent.
(or, There is no relationship between A and B.)
HA : Variables A and B are not independent.
(or, There is a relationship between A and B.)
To carry out the test, therefore, you compare the observed data to the
counts you'd expect to see if the distributions were completely unrelated
(satisfying the assumption that H0 is true). An overall large departure from
the expected frequencies will result in a small p-value and thus provide evidence against the null.
So, how are such data best presented? For two categorical variables, a
two-dimensional structure is appropriate; in R, this is a standard matrix.
For example, suppose some dermatologists at a certain clinical practice
are interested in their successes in treating a common skin affliction. Their
records show N = 355 patients have been treated at their clinic using one of
four possible treatments-a course of tablets, a series of injections, a laser
treatment, and an herbal-based remedy. The level of success in curing the
affliction is also recorded-none, partial success, and full success. The data
are given in the constructed matrix skin.
R> skin <- matrix(c(20,32,8,52,9,72,8,32,16,64,30,12),4,3,
                  dimnames=list(c("Injection","Tablet","Laser","Herbal"),
                                c("None","Partial","Full")))
Hypothesis Testing 415
R> skin
None Partial Full
Injection 20 9 16
Tablet 32 72 64
Laser 8 8 30
Herbal 52 32 12
A two-dimensional table presenting frequencies in this fashion is called a
contingency table.
Calculation: Chi-Squared Test of Independence
To compute the test statistic, presume data are presented as a kr × kc contingency table, in other words, a matrix of counts, based on two categorical
variables (both mutually exclusive and exhaustive). The focus of the test is
the way in which the frequencies of N observations between the kr levels
of the "row" variable and the kc levels of the "column" variable are jointly
distributed. The test statistic ??
2
is given with
??
2 =
  X
kr
i=1
X
kc
j=1
(O[i, j] ??? E[i, j])
2
E[i, j]
, (18.12)
where O[i, j] is the observed count and E[i, j] is the expected count at row
position i and column position j. Each E[i, j] is found as the sum total of row
i multiplied by the sum total of column j, all divided by N.
E[i, j] =
  Pkr
u=1
O[u, j]

×
Pkc
v=1
O[i,v]

N
(18.13)
The result, ??
2
, follows a chi-squared distribution with ?? = (kr ??? 1) ×
(kc ??? 1) degrees of freedom. Again, the p-value is always an upper-tailed
area, and you can consider the test valid with the satisfaction of the condition that at least 80 percent of the E[i, j] are at least 5.
For this calculation, it's important to note the following:
  . It's not necessary to assume that kr = kc .
. The functionality of Equation (18.12) is the same as that of (18.11)-an
overall sum involving the squared differences between the observed and
expected values of each cell.
. The double-sum in (18.12) just represents the total sum over all the
cells, in the sense that you can compute the total sample size N with
Pkr
i=1
Pkc
j=1
O[i, j].
. A rejected H0 doesn't tell you about the nature of how the frequencies
depend on one another, just that there is evidence to suggest that some
kind of dependency between the two categorical variables exists.
416 Chapter 18
Continuing with the example, the dermatologists want to determine
whether their records suggest there is statistical evidence to indicate some
relationship between the type of treatment and the level of success in curing
the skin affliction. For convenience, store the total number of categories kr
and kc for the row and column variables, respectively.
R> kr <- nrow(skin)
R> kc <- ncol(skin)
You have the O[i, j] in skin, so now you must now compute the E[i, j]. In
light of Equation (18.13), which deals with row and column sums, you can
evaluate these using the built-in rowSums and colSums functions.
R> rowSums(skin)
Injection Tablet Laser Herbal
45 168 46 96
R> colSums(skin)
None Partial Full
112 121 122
These results indicate the totals in each group, regardless of the
other variable. To get the expected counts for all cells of the matrix, Equation (18.13) requires each row sum to be multiplied by each column sum
once. You could write a for loop, but this would be inefficient and rather
inelegant. It is better to use rep with the optional each argument (refer to
                                                                    Section 2.3.2). By repeating each element of the column totals (level of success) four times, you can then use vector-oriented behavior to multiply that
repeated vector by the shorter vector produced by rowSums. You can then call
sum(skin) to divide this by N and rearrange it into a matrix. The following
lines show how this example works step-by-step:
  R> rep(colSums(skin),each=kr)
None None None None Partial Partial Partial Partial Full
112 112 112 112 121 121 121 121 122
Full Full Full
122 122 122
R> rep(colSums(skin),each=kr)*rowSums(skin)
None None None None Partial Partial Partial Partial Full
5040 18816 5152 10752 5445 20328 5566 11616 5490
Full Full Full
20496 5612 11712
R> rep(colSums(skin),each=kr)*rowSums(skin)/sum(skin)
None None None None Partial Partial Partial Partial
14.19718 53.00282 14.51268 30.28732 15.33803 57.26197 15.67887 32.72113
Full Full Full Full
15.46479 57.73521 15.80845 32.99155
R> skin.expected <- matrix(rep(colSums(skin),each=kr)*rowSums(skin)/sum(skin),
                           nrow=kr,ncol=kc,dimnames=dimnames(skin))
Hypothesis Testing 417
R> skin.expected
None Partial Full
Injection 14.19718 15.33803 15.46479
Tablet 53.00282 57.26197 57.73521
Laser 14.51268 15.67887 15.80845
Herbal 30.28732 32.72113 32.99155
Note that all the expected values are greater than 5, as preferred.
It's best to construct a single object to hold the results of the different
stages of calculations leading to the test statistic, as you did for the onedimensional example. Since each stage is a matrix, you can bind the relevant matrices together with cbind and produce an array of the appropriate
dimensions (refer to Section 3.4 for a refresher).
R> skin.array <- array(data=cbind(skin,skin.expected,
                                  (skin-skin.expected)^2/skin.expected),
                       dim=c(kr,kc,3),
                       dimnames=list(dimnames(skin)[[1]],dimnames(skin)[[2]],
                                     c("O[i,j]","E[i,j]",
                                       "(O[i,j]-E[i,j])^2/E[i,j]")))
R> skin.array
, , O[i,j]
None Partial Full
Injection 20 9 16
Tablet 32 72 64
Laser 8 8 30
Herbal 52 32 12
, , E[i,j]
None Partial Full
Injection 14.19718 15.33803 15.46479
Tablet 53.00282 57.26197 57.73521
Laser 14.51268 15.67887 15.80845
Herbal 30.28732 32.72113 32.99155
, , (O[i,j]-E[i,j])^2/E[i,j]
None Partial Full
Injection 2.371786 2.6190199 0.01852279
Tablet 8.322545 3.7932587 0.67978582
Laser 2.922614 3.7607992 12.74002590
Herbal 15.565598 0.0158926 13.35630339
418 Chapter 18
The final steps are easy-the test statistic given by (18.12) is just the
grand total of all elements of the matrix that is the third layer of skin.array.
R> X2 <- sum(skin.array[,,3])
R> X2
[1] 66.16615
The corresponding p-value for this test of independence is as follows:
  R> 1-pchisq(X2,df=(kr-1)*(kc-1))
[1] 2.492451e-12
Recall that the relevant degrees of freedom are defined as ?? = (kr ??? 1) ×
(kc ??? 1).
The extremely small p-value provides strong evidence against the null
hypothesis. The appropriate conclusion would be to reject H0 and state that
there does appear to be a relationship between the type of treatment for the
skin affliction and the level of success in curing it.
R Function: chisq.test
Yet once more, no section in this chapter would be complete without showcasing the built-in functionality R possesses for these fundamental procedures. The default behavior of chisq.test, when supplied a matrix as x, is
to perform a chi-squared test of independence with respect to the row and
column frequencies-just as performed manually here for the skin affliction
example. The following result easily confirms your previous calculations:
  R> chisq.test(x=skin)
Pearson's Chi-squared test
data: skin
X-squared = 66.1662, df = 6, p-value = 2.492e-12
Exercise 18.4
HairEyeColor is a ready-to-use data set in R that you haven't yet come
across. This 4 × 4 × 2 array provides frequencies of hair and eye colors
of 592 statistics students, split by sex (Snee, 1974).
a. Perform and interpret, at a significance level of ?? = 0.01, a chisquared test of independence for hair against eye color for all
students, regardless of their sex.
Hypothesis Testing 419
In Exercise 8.1 on page 161, you accessed the Duncan data set of the
contributed package car, which contains markers of job prestige
collected in 1950. Install the package if you haven't already and load
the data frame.
b. The first column of Duncan is the variable type, recording the
type of job as a factor with three levels: prof (professional or
managerial), bc (blue collar), and wc (white collar). Construct
appropriate hypotheses and perform a chi-squared GOF test to
determine whether the three job types are equally represented in
the data set.
i. Interpret the resulting p-value with respect to a significance
level of ?? = 0.05.
ii. What conclusion would you reach if you used a significance
level of ?? = 0.01?
18.5 Errors and Power
In discussing all these forms of statistical hypothesis testing, there has been
one common thread: the interpretation of a p-value and what it tells you
about your problem in terms of the hypotheses. Frequentist statistical
hypothesis testing is ubiquitous in many fields of research, so it is important to at least briefly explore directly related concepts.
18.5.1 Hypothesis Test Errors
Hypothesis testing is performed with the objective of obtaining a p-value in
order to quantify evidence against the null statement H0. This is rejected in
favor of the alternative, HA, if the p-value is itself less than a predefined significance level ??, which is conventionally 0.05 or 0.01. As touched upon, this
approach is justifiably criticized since the choice of ?? is essentially arbitrary;
a decision to reject or retain H0 can change depending solely upon the ??
value.
Consider for the moment, given a specific test, what the correct outcome
is. If H0 is really true, then you'd want to retain it. If HA is really true, you'd
want to reject the null. This "truth," one way or another, is impossible to
know in practice. That being said, it's useful to consider in a theoretical
sense just how good (or bad) a given hypothesis test is at yielding a result
that leads to the correct conclusion.
To be able to test the validity of your rejection or retention of the null
hypothesis, you must be able to identify two kinds of errors:
. A Type I error occurs when you incorrectly reject a true H0. In any given
hypothesis test, the probability of a Type I error is equivalent to the significance level ??.
420 Chapter 18
. A Type II error occurs when you incorrectly retain a false H0 (in other
words, fail to accept a true HA). Since this depends upon what the true
HA actually is, the probability of committing such an error, labeled ??, is
not usually known in practice.
18.5.2 Type I Errors
If your p-value is less than ??, you reject the null statement. If the null is
really true, though, the ?? directly defines the probability that you incorrectly
reject it. This is referred to as a Type I error.
Figure 18-2 provides a conceptual illustration of a Type I error probability for a supposed hypothesis test of a sample mean, where the hypotheses
are set up as H0 : µ = µ0 and HA : µ > µ0.
Figure 18-2: A conceptual diagram of the Type I error
probability ??
The null hypothesis distribution is centered on the null value µ0; the
alternative hypothesis distribution is centered to its right at some mean µA
in Figure 18-2. As you can see, if the null hypothesis is really true, then the
probability it is incorrectly rejected for this test will be equal to the significance level ??, located in the upper tail of the null distribution.
Simulating Type I Errors
To demonstrate the Type I error rate via numerical simulation (here, this
refers to randomly generating hypothetical data samples), you can write
code that does the equivalent of repeating a hypothesis test under known
conditions. So that you can use this code multiple times, in the R script editor define the following function:
typeI.tester <- function(mu0,sigma,n,alpha,ITERATIONS=10000){
pvals <- rep(NA,ITERATIONS)
for(i in 1:ITERATIONS){
temporary.sample <- rnorm(n=n,mean=mu0,sd=sigma)
Hypothesis Testing 421
temporary.mean <- mean(temporary.sample)
temporary.sd <- sd(temporary.sample)
pvals[i] <- 1-pt((temporary.mean-mu0)/(temporary.sd/sqrt(n)),df=n-1)
}
return(mean(pvals<alpha))
}
The typeI.tester function is designed to generate ITERATIONS samples
from a particular normal distribution. With each sample, you'll perform an
upper-tailed test of the mean (refer to Section 18.2.1) in the spirit of Figure 18-2, assuming the hypotheses of H0 : µ = µ0 and HA : µ > µ0.
You can decrease ITERATIONS to generate fewer entire samples, and this
will speed up computation time but will result in simulated rates that are
more variable. Each entire sample of size n of hypothetical raw measurements is generated using rnorm with the mean equal to the mu0 argument
(and standard deviation equal to the sigma argument). The desired significance level is set by alpha. In the for loop, the sample mean and sample
standard deviation are calculated for each generated sample.
Were each sample subjected to a "real" hypothesis test, the p-value
would be taken from the right-hand area of the t-distribution with n-1
degrees of freedom (using pt), with respect to the standardized test statistic given earlier in Equation (18.2).
The calculated p-value, at each iteration, is stored in a predefined vector pvals. The logical vector pvals<alpha therefore contains corresponding
TRUE/FALSE values; the former logical value flags rejection of the null hypothesis, and the latter flags retention. The Type I error rate is determined by
calling mean on that logical vector, which yields the proportion of TRUEs (in
other words, the overall proportion of "null hypothesis rejections") arising
from the simulated samples. Remember, the samples are generated randomly, so your results are liable to change slightly each time you run the
function.
This function works because, by definition of the problem, the samples
that are being generated come from a distribution that truly has the mean
set at the null value, in other words, µA = µ0. Therefore, any statistical rejection of this statement, obtained with a p-value less than the significance level
??, is clearly incorrect and is purely a result of random variation.
To try this, import the function and execute it generating the default
ITERATIONS=10000 samples. Use the standard normal as the null (and "true" in
this case!) distribution; make each sample of size 40 and set the significance
level at the conventional ?? = 0.05. Here's an example:
R> typeI.tester(mu0=0,sigma=1,n=40,alpha=0.05)
[1] 0.0489
422 Chapter 18
This indicates that 10,000 × 0.0489 = 489 of the samples taken yielded a
corresponding test statistic that provided a p-value, which would incorrectly
result in rejection of H0. This simulated Type I error rate lies close to the
preset alpha=0.05.
Here's another example, this time for nonstandard normal data samples
with ?? = 0.01:
R> typeI.tester(mu0=-4,sigma=0.3,n=60,alpha=0.01)
[1] 0.0108
Note that again, the numerically simulated rate of Type I error reflects
the significance level.
These results are not difficult to understand theoretically-if the true
distribution does indeed have a mean equal to the null value, you'll naturally
observe those "extreme" test statistic values in practice at a rate equal to ??.
The catch, of course, is that in practice the true distribution is unknown,
highlighting once more the fact that a rejection of any H0 can never be
interpreted as proof of the truth of HA. It might simply be that the sample
you observed followed the null hypothesis but produced an extreme test
statistic value by chance, however small that chance might be.
Bonferroni Correction
The fact that Type I errors naturally occur because of random variation is
particularly important and leads us to consider the multiple testing problem. If
you're conducting many hypothesis tests, you should be cautious in simply
reporting the "number of statistically significant outcomes"-as you increase
the number of hypothesis tests, you increase the chance of receiving an erroneous result. In, say, 20 tests conducted under ?? = 0.05, on average one will
be a so-called false positive; if you conduct 40 or 60 tests, you are inevitably
more likely to find more false positives.
When several hypothesis tests are conducted, you can curb the multiple
testing problem with respect to committing a Type I error by using the Bonferroni correction. The Bonferroni correction suggests that when performing
a total of N independent hypothesis tests, each under a significance level of
??, you should instead use ??B = ??/N for any interpretation of statistical significance. Be aware, however, that this correction to the level of significance
represents the simplest solution to the multiple testing problem and can be
criticized for its conservative nature, which is potentially problematic when
N is large.
The Bonferroni and other corrective measures were developed in an
attempt to formalize remedies to making a Type I error in multiple tests. In
general, though, it suffices to be aware of the possibility that H0 may be true,
even if the p-value is considered small.
Hypothesis Testing 423
18.5.3 Type II Errors
The issues with Type I errors might suggest that it's desirable to perform a
hypothesis test with a smaller ?? value. Unfortunately, it's not quite so simple;
reducing the significance level for any given test leads directly to an increase
in the chance of committing a Type II error.
A Type II error refers to incorrect retention of the null hypothesis-in
other words, obtaining a p-value greater than the significance level when it's
the alternative hypothesis that's actually true. For the same scenario you've
been looking at so far (an upper-tailed test for a single sample mean), Figure 18-3 illustrates the probability of a Type II error, shaded and denoted ??.
Figure 18-3: A conceptual diagram of the Type II error
probability ??
It's not as easy to find ?? as it is to find the probability of making a Type I
error because ?? depends, among other things, on what the true value of
µA is (which in general you won't know). If µA is closer to the hypothesized
null value of µ0, you can imagine the alternative distribution in Figure 18-3
translating (shifting) to the left, resulting in an increase in ??. Similarly, staying with Figure 18-3, imagine decreasing the significance level ??. Doing so
means the vertical dashed line (denoting the corresponding critical value)
moves to the right, also increasing the shaded area of ??. Intuitively, this
makes sense-the closer the true alternative value is to the null and/or the
smaller the significance level, the harder HA is to detect by rejection of H0.
As noted, ?? usually can't be calculated in practice because of the need
to know what the true distribution actually is. This quantity is, however, useful in giving you an idea of how prone a test is to the incorrect retention of
a null hypothesis under particular conditions. Suppose, for example, you're
performing a one-sample t-test for H0 : µ = µ0 and HA : µ > µ0 with µ0 = 0
but that the (true) alternative distribution of the raw measurements has
mean µA = 0.5 and standard deviation ?? = 1. Given a random sample of
size n = 30 and using ?? = 0.05, what is the probability of committing a
Type II error in any given hypothesis test (using the same standard deviation
424 Chapter 18
for the null distribution)? To answer this, look again at Figure 18-3; you
need the critical value marked off by the significance level (the dashed
vertical line). If you assume ?? is known, then the sampling distribution of
interest will be normal with mean µ0 = 0 and a standard error of 1/
???
30 (see
Section 17.1.1). Therefore, with an upper-tail area of 0.05, you can find the
critical value with the following:
R> critval <- qnorm(1-0.05,mean=0,sd=1/sqrt(30))
R> critval
[1] 0.3003078
This represents the vertical dashed line in this specific setting (see
Section 16.2.2 for a refresher on use of qnorm). The Type II error in this
example is found as the left-hand tail area under the alternative, "true" distribution, from that critical value:
R> pnorm(critval,mean=0.5,sd=1/sqrt(30))
[1] 0.1370303
From this, you can see that a hypothesis test under these conditions has
roughly a 13.7 percent chance of incorrect retention of the null.
Simulating Type II Errors
Simulation is especially useful here. In the editor, consider the function
typeII.tester defined as follows:
typeII.tester <- function(mu0,muA,sigma,n,alpha,ITERATIONS=10000){
pvals <- rep(NA,ITERATIONS)
for(i in 1:ITERATIONS){
temporary.sample <- rnorm(n=n,mean=muA,sd=sigma)
temporary.mean <- mean(temporary.sample)
temporary.sd <- sd(temporary.sample)
pvals[i] <- 1-pt((temporary.mean-mu0)/(temporary.sd/sqrt(n)),df=n-1)
}
return(mean(pvals>=alpha))
}
This function is similar to typeI.tester. The null value, standard deviation of raw measurements, sample size, significance level, and number of
iterations are all as before. Additionally, you now have muA, providing the
"true" mean µA under which to generate the samples. Again, at each iteration, a random sample of size n is generated, its mean and standard deviation are calculated, and the appropriate p-value for the test is computed
using pt from the usual standardized test statistic with df=n-1. (Remember,
since you're estimating the true standard deviation of the measurements
?? with the sample standard deviation s, it's technically correct to use the
Hypothesis Testing 425
t-distribution.) Following completion of the for loop, the proportion of
p-values that were greater than or equal to the significance level alpha is
returned.
After importing the function into the workspace, you can simulate ?? for
this test.
R> typeII.tester(mu0=0,muA=0.5,sigma=1,n=30,alpha=0.05)
[1] 0.1471
My result indicates something close to the theoretical ?? evaluated previously, albeit slightly larger because of the additional uncertainty that is
naturally present when using a t-based sampling distribution instead of a
normal. Again, each time you run typeII.tester, the results will vary slightly
since everything is based on randomly generated hypothetical data samples.
Turning your attention to Figure 18-3, you can see (in line with a comment made earlier) that if, in an effort to decrease the chance of a Type I
error, you use ?? = 0.01 instead of 0.05, the vertical line moves to the right,
thereby increasing the probability of a Type II error, with all other conditions being held constant.
R> typeII.tester(mu0=0,muA=0.5,sigma=1,n=30,alpha=0.01)
[1] 0.3891
Other Influences on the Type II Error Rate
The significance level isn't the only contributing factor in driving ??. Keeping ?? at 0.01, this time see what happens if the standard deviation of the raw
measurements is increased from ?? = 1 to ?? = 1.1 and then ?? = 1.2.
R> typeII.tester(mu0=0,muA=0.5,sigma=1.1,n=30,alpha=0.01)
[1] 0.4815
R> typeII.tester(mu0=0,muA=0.5,sigma=1.2,n=30,alpha=0.01)
[1] 0.5501
Increasing the variability of the measurements, without touching anything else in the scenario, also increases the chance of a Type II error. You
can imagine the curves in Figure 18-3 becoming flatter and more widely dispersed owing to a larger standard error of the mean, which would result in
more probability weight in the left-hand tail marked off by the critical value.
Conversely, if the variability of the raw measurements is smaller, then the
sampling distributions of the sample mean will be taller and skinnier, meaning a reduction in ??.
A smaller or larger sample size will have a similar impact. Located in
the denominator of the standard error formula, a smaller n will result in a
larger standard error and hence that flatter curve and an increased ??;
a larger sample size will have the opposite effect. If you remain with the
latest values of µ0 = 0, µA = 0.5, ?? = 1.2, and ?? = 0.01, note that reducing
426 Chapter 18
the sample size to 20 (from 30) results in an increased simulated Type II
error rate compared with the most recent result of 0.5501, but increasing
the sample size to 40 improves the rate.
R> typeII.tester(mu0=0,muA=0.5,sigma=1.2,n=20,alpha=0.01)
[1] 0.7319
R> typeII.tester(mu0=0,muA=0.5,sigma=1.2,n=40,alpha=0.01)
[1] 0.4219
Finally, as noted at the beginning of the discussion, the specific value of
µA itself affects ?? just as you'd expect. Again, keeping the latest values for
all other components, which resulted in my case in ?? = 0.4219, note that
shifting the "true" mean closer to µ0 by changing from µA = 0.5 to µA = 0.4
means the probability of committing a Type II error is increased; the opposite is true if the difference is increased to µA = 0.6.
R> typeII.tester(mu0=0,muA=0.4,sigma=1.2,n=40,alpha=0.01)
[1] 0.6147
R> typeII.tester(mu0=0,muA=0.6,sigma=1.2,n=40,alpha=0.01)
[1] 0.2287
To summarize, although these simulated rates have been applied to
the specific situation in which the hypothesis test is an upper-tailed test for
a single mean, the general concepts and ideas discussed here hold for any
hypothesis test. It's easy to establish that the Type I error rate matches the
predefined significance level and so can be decreased by reducing ??. In
contrast, controlling the Type II error rate is a complex balancing act that
can involve sample size, significance level, observation variability, and magnitude of the difference between the true value and the null. This problem is
largely academic since the "truth" is typically unknown in practice. However,
the Type II error rate's direct relationship to statistical power often plays a
critical role in preparing for data collection, especially when you're considering sample size requirements, as you'll see in the next section.
Exercise 18.5
a. Write a new version of typeI.tester called typeI.mean. The new
function should be able to simulate the Type I error rate for
tests of a single mean in any direction (in other words, oneor two-sided). The new function should take an additional
argument, test, which takes a character string "less", "greater",
or "two.sided" depending on the type of desired test. You can
achieve this by modifying typeI.tester as follows:
- Instead of calculating and storing the p-values directly in the
for loop, simply store the test statistic.
Hypothesis Testing 427
- When the loop is complete, set up stacked if-else statements
that cater to each of the three types of test, calculating the
p-value as appropriate.
- For the two-sided test, remember that the p-value is defined
as twice the area "more extreme" than the null. Computationally, this means you must use the upper-tail area if the
test statistic is positive and the lower-tail area otherwise.
If this area is less than half of ?? (since it is subsequently
multiplied by 2 in a "real" hypothesis test), then a rejection
of the null should be flagged.
- If the value of test is not one of the three possibilities, the
function should throw an appropriate error using stop.
i. Experiment with your function using the first example setting in the text with µ0 = 0, ?? = 1, n = 40, and ?? = 0.05.
Call typeI.mean three times, using each of the three possible
options for test. You should find that all simulated results sit
close to 0.05.
ii. Repeat (i) using the second example setting in the text with
µ0 = ???4, ?? = 0.3, n = 60, and ?? = 0.01. Again, you should
find that all simulated results sit close to the value of ??.
b. Modify typeII.tester in the same way as you did typeI.tester; call
the new function typeII.mean. Simulate the Type II error rates
for the following hypothesis tests. As per the text, assume µA,
??, ??, and n denote the true mean, standard deviation of raw
observations, significance level, and sample size, respectively.
i. H0 : µ = ???3.2; HA : µ , ???3.2
with µA = ???3.3, ?? = 0.1, ?? = 0.05, and n = 25.
ii. H0 : µ = 8994; HA : µ < 8994
with µA = 5600, ?? = 3888, ?? = 0.01, and n = 9.
iii. H0 : µ = 0.44; HA : µ > 0.44
with µA = 0.4, ?? = 2.4, ?? = 0.05, and n = 68.
18.5.4 Statistical Power
For any hypothesis test, it is useful to consider its potential statistical power.
Power is the probability of correctly rejecting a null hypothesis that is untrue.
For a test that has a Type II error rate of ??, the statistical power is found
simply with 1 ??? ??. It's desirable for a test to have a power that's as high as
possible. The simple relationship with the Type II error probability means
that all factors impacting the value of ?? also directly affect power.
For the same one-sided H0 : µ = µ0 and HA : µ > µ0 example discussed
in the previous section, Figure 18-4 shades the power of the test-the complement to the Type II error rate. By convention, a hypothesis test that has a
power greater than 0.8 is considered statistically powerful.
428 Chapter 18
Figure 18-4: A conceptual diagram of statistical power 1 ??? ??
You can numerically evaluate power under specific testing conditions
using simulation. For the previous discussion on Type II errors, you're able
to subtract all simulated results of ?? from 1 to evaluate the power of that particular test. For example, the power of detection of µA = 0.6 when µ0 = 0,
taking samples of size n = 40 and with ?? = 1.2 and ?? = 0.01, is simulated as
1 ??? 0.2287 = 0.7713 (using my most recent result of ?? from earlier). This
means there's roughly a 77 percent chance of correctly detecting the true
mean of 0.6 in a hypothesis test based on a sample of measurements generated under those conditions.
Researchers are often interested in the relationship between power and
sample size (though it is important to bear in mind that this is only one of
the influential ingredients in the determination of power). Before you begin
to collect data to examine a particular hypothesis, you might have an idea of
the potential true value of the parameter of interest from past research or
pilot studies. This is useful in helping to determine your sample size, such as
in helping to answer questions like "How big does my sample need to be in
order to be able to conduct a statistically powerful test to correctly reject H0,
if the true mean is actually µA?"
Simulating Power
For the most recent testing conditions, with a sample size of n = 40, you've
seen that there's a power of around 0.77 of detecting µA = 0.6. For the
purposes of this example, let's say you want to find how much you should
increase n by in order to conduct a statistically powerful test. To answer this,
define the following function power.tester in the editor:
power.tester <- function(nvec,...){
nlen <- length(nvec)
result <- rep(NA,nlen)
pbar <- txtProgressBar(min=0,max=nlen,style=3)
for(i in 1:nlen){
Hypothesis Testing 429
result[i] <- 1-typeII.tester(n=nvec[i],...)
setTxtProgressBar(pbar,i)
}
close(pbar)
return(result)
}
The power.tester function uses the typeII.tester function defined in Section 18.5.3 to evaluate the power of a given upper-tailed hypothesis test of
a single sample mean. It takes a vector of sample sizes supplied as the nvec
argument (you pass all other arguments to typeII.tester using an ellipsis-
refer to Section 11.2.4). A for loop defined in power.tester cycles through
the entries of nvec one at a time, simulates the power for each sample size,
and stores them in a corresponding vector that's returned to the user.
Remember, through typeII.tester, this function is using random generation of hypothetical data samples, so there may be some fluctuation in the
results you observe each time you run power.tester.
There can be a slight delay when evaluating the power for many individual sample sizes, so this function also provides a good opportunity to showcase a progress bar in a practical implementation (refer to Section 12.2.1 for
details).
Set up the following vector, which uses the colon operator (see Section 2.3.2) to construct a sequence of integers between 5 and 100 inclusive
for the sample sizes to be examined:
R> sample.sizes <- 5:100
Importing the power.tester function, you can then simulate the power
for each of these integers for this particular test (ITERATIONS is halved to 5000
to reduce the overall completion time).
R> pow <- power.tester(nvec=sample.sizes,
mu0=0,muA=0.6,sigma=1.2,alpha=0.01,ITERATIONS=5000)
|====================================================================| 100%
R> pow
[1] 0.0630 0.0752 0.1018 0.1226 0.1432 0.1588 0.1834 0.2162 0.2440 0.2638
[11] 0.2904 0.3122 0.3278 0.3504 0.3664 0.3976 0.4232 0.4478 0.4680 0.4920
[21] 0.5258 0.5452 0.5552 0.5616 0.5916 0.6174 0.6326 0.6438 0.6638 0.6844
[31] 0.6910 0.7058 0.7288 0.7412 0.7552 0.7718 0.7792 0.7950 0.8050 0.8078
[41] 0.8148 0.8316 0.8480 0.8524 0.8600 0.8702 0.8724 0.8800 0.8968 0.8942
[51] 0.8976 0.9086 0.9116 0.9234 0.9188 0.9288 0.9320 0.9378 0.9370 0.9448
[61] 0.9436 0.9510 0.9534 0.9580 0.9552 0.9648 0.9656 0.9658 0.9684 0.9756
[71] 0.9742 0.9770 0.9774 0.9804 0.9806 0.9804 0.9806 0.9854 0.9848 0.9844
[81] 0.9864 0.9886 0.9890 0.9884 0.9910 0.9894 0.9906 0.9930 0.9926 0.9938
[91] 0.9930 0.9946 0.9948 0.9942 0.9942 0.9956
As expected, the power of detection rises steadily as n increases; the
conventional cutoff of 80 percent is visible in these results as lying between
430 Chapter 18
0.7950 and 0.8050. If you don't want to identify the value visually, you can
find which entry of sample.sizes corresponds to the 80 percent cutoff by first
using which to identify the indexes of pow that are at least 0.8 and then returning the lowest value in that category with min. The identified index may then
be specified in square brackets to sample.sizes to give you the value of n that
corresponds to that simulated power (0.8050 in this case). These commands
can be nested as follows:
R> minimum.n <- sample.sizes[min(which(pow>=0.8))]
R> minimum.n
[1] 43
The result indicates that if your sample size is at least 43, a hypothesis
test under these particular conditions should be statistically powerful (based
on the randomly simulated output in pow in this instance).
What if the significance level for this test were relaxed? Say you wanted
to conduct the test (still upper-tailed under the condition of µ0 = 0,
µA = 0.6, and ?? = 1.2) using a significance level of ?? = 0.05 rather than
0.01. If you look again at Figure 18-4, this alteration means the vertical line
(critical value) moves to the left, decreasing ?? and so increasing power. That
would therefore suggest you'd require a smaller sample size than earlier, in
other words, n < 43, in order to perform a statistically powerful test when ??
is increased.
To simulate this situation for the same range of sample sizes and store
the resulting powers in pow2, examine the following:
R> pow2 <- power.tester(nvec=sample.sizes,
mu0=0,muA=0.6,sigma=1.2,alpha=0.05,ITERATIONS=5000)
|====================================================================| 100%
R> minimum.n2 <- sample.sizes[min(which(pow2>0.8))]
R> minimum.n2
[1] 27
This result indicates a sample size of at least 27 is required, which is a
noticeable reduction from the 43 noted if ?? = 0.01. However, relaxing ??
means an increased risk of committing a Type I error!
Power Curves
For comparison, you can plot your simulated powers as a kind of power
curve using both pow and pow2 with the following code:
R> plot(sample.sizes,pow,xlab="sample size n",ylab="simulated power")
R> points(sample.sizes,pow2,col="gray")
R> abline(h=0.8,lty=2)
R> abline(v=c(minimum.n,minimum.n2),lty=3,col=c("black","gray"))
R> legend("bottomright",legend=c("alpha=0.01","alpha=0.05"),
col=c("black","gray"),pch=1)
Hypothesis Testing 431
My particular image is given in Figure 18-5. A horizontal line marks off
the power of 0.8, and the vertical line marks the minimum sample size values
identified and stored in minimum.n and minimum.n2. As a final touch, a legend is
added to reference the ?? values of each curve.
Figure 18-5: Simulated power curves for the upper-tailed
hypothesis test of a single sample mean
The curves themselves indicate exactly what you'd expect-the power
of detection increases as the sample size is incremented. You can also note
the flattening off occurring as the power rises ever closer to the "perfect"
rate of 1, which is typical of a power curve. For ?? = 0.05, the curve sits almost
consistently above the curve for ?? = 0.01, though the difference looks negligible as n rises above 75 or so.
The preceding discussion of errors and power highlights the need for
care in interpreting the results of even the most basic of statistical tests.
A p-value is merely a probability, and as such, no matter how small it may
be in any circumstance, it can never prove or disprove a claim on its own.
Issues surrounding the quality of a hypothesis test (parametric or otherwise)
should be considered, though this is arguably difficult in practice. Nevertheless, an awareness of Type I and Type II errors, as well as the concept of
statistical power, is extremely useful in the implementation and appraisal of
any formal statistical testing procedure.
432 Chapter 18
Exercise 18.6
a. For this exercise you'll need to have written typeII.mean from
Exercise 18.5 (b). Using this function, modify power.tester so
that a new function, power.mean, calls typeII.mean instead of calling
typeII.tester.
i. Confirm that the power of the test given by H0 : µ = 10;
HA : µ , 10, with µA = 10.5, ?? = 0.9, ?? = 0.01, and n = 50, is
roughly 88 percent.
ii. Remember the hypothesis test in Section 18.2.1 for the mean
net weight of an 80-gram pack of snacks, based on the n = 44
observations provided in the snack vector. The hypotheses
were as follows:
H0 : µ = 80
HA : µ < 80
If the true mean is µA = 78.5 g and the true standard
deviation of the weights is ?? = 3.1 g, use power.mean to determine whether the test is statistically powerful, assuming
?? = 0.05. Does your answer to this change if ?? = 0.01?
b. Staying with the snacks hypothesis test, using the sample.sizes vector from the text, determine the minimum sample size required
for a statistically powerful test using both ?? = 0.05 and ?? = 0.01.
Produce a plot showing the two power curves.
Important Code in This Chapter
Function/operator Brief description First occurrence
t.test One- and two-sample t-test Section 18.2.1, p. 391
prop.test One- and two-sample Z-test Section 18.3.1, p. 405
pchisq ??
2
cumulative problems Section 18.4.1, p. 414
chisq.test ??
2
test of distribution/independence Section 18.4.1, p. 414
rowSums Matrix row totals Section 18.4.2, p. 417
colSums Matrix column totals Section 18.4.2, p. 417
Hypothesis Testing 433
19
ANALYS IS OF VAR IANCE
Analysis of variance (ANOVA), in its simplest
form, is used to compare multiple means
in a test for equivalence. In that sense, it's a
straightforward extension of the hypothesis test
comparing two means. There's a continuous variable
from which the means of interest are calculated, and
there's at least one categorical variable that tells you how to define the
groups for those means. In this chapter, you'll explore the ideas surrounding ANOVA and look at comparing means first split by one categorical
variable (one-way ANOVA) and then split by multiple categorical variables
(multiple-factor ANOVA).
19.1 One-Way ANOVA
The simplest version of ANOVA is referred to as one-way or one-factor analysis.
Simply put, the one-way ANOVA is used to test two or more means for equality. Those means are split by a categorical group or factor variable. ANOVA is
often used to analyze experimental data to assess the impact of an intervention. You might, for example, be interested in comparing the mean weights
of the chicks in the built-in chickwts data set, split according to the different
food types they were fed.
19.1.1 Hypotheses and Diagnostic Checking
Say you have a categorical-nominal variable that splits a total of N numeric
observations into k distinct groups, where k ??? 2. You're looking to statistically compare the k groups' means, µ1,. . . , µk , to see whether they can be
claimed to be equal. The standard hypotheses are as follows:
H0 : µ1 = µ2 = . . . = µk
HA : µ1, µ2, . . . , µk are not all equal
(alternatively, at least one mean differs).
In fact, when k = 2, the two-sample t-test is equivalent to ANOVA; for
that reason, ANOVA is most frequently employed when k > 2.
The following assumptions need to be satisfied in order for the results of
the basic one-way ANOVA test to be considered reliable:
Independence The samples making up the k groups must be independent of one another, and the observations in each group must be
independent and identically distributed (iid).
Normality The observations in each group should be normally distributed, or at least approximately so.
Equality of variances The variance of the observations in each group
should be equal, or at least approximately so.
If the assumptions of equality of variances or normality are violated, it
doesn't necessarily mean your results will be completely worthless, but it will
impact the overall effectiveness of detecting a true difference in the means
(refer to the discussion on statistical power in Section 18.5.4). It's always a
good idea to assess the validity of these assumptions before using ANOVA;
I'll do this informally for the upcoming example.
It's also worth noting that you don't need to have an equal number of
observations in each group to perform the test (in which case it is referred
to as unbalanced). However, having unbalanced groups does render the test
more sensitive to potentially detrimental effects if your assumptions of equality of variances and normality are not sound.
Let's return to the chickwts data for the example-the weights of chicks
based on k = 6 different feeds. You're interested in comparing the mean
weights according to feed type to see whether they're all equal. Use table
to summarize the six sample sizes and use tapply (see, for example, Section 13.2.1) to get each group mean, as follows:
R> table(chickwts$feed)
casein horsebean linseed meatmeal soybean sunflower
12 10 12 11 14 12
R> chick.means <- tapply(chickwts$weight,INDEX=chickwts$feed,FUN=mean)
R> chick.means
casein horsebean linseed meatmeal soybean sunflower
323.5833 160.2000 218.7500 276.9091 246.4286 328.9167
436 Chapter 19
Your skills from Section 14.3.2 allow you to produce side-by-side boxplots of the distributions of weights. The next two lines give you the plot on
the left of Figure 19-1:
R> boxplot(chickwts$weight~chickwts$feed)
R> points(1:6,chick.means,pch=4,cex=1.5)
Figure 19-1: Exploring the chickwts data. Left: Side-by-side boxplots of chick weight split
by feed type, with the mean marked by ×. Right: Normal QQ plot of the mean-centered
data of each feed group.
Because boxplots display the median, not the mean, the second line of
code adds the feed-specific means (stored in the chick.means object you just
created) to each box using points.
Inspecting the left plot of Figure 19-1, it certainly looks as though
there's a difference in the mean weights. Is any apparent difference statistically significant, though? To find out, the ANOVA test for this example
concerns the following hypotheses:
H0 : µcasein = µhorsebean = µlinseed = µmeatmeal = µsoybean = µsunflower
HA : The means are not all equal.
Assuming independence of the data, before implementing the test, you
must first check that the other assumptions are valid. To examine equality of
variances, you can use the same informal rule of thumb as used in the twosample t-test. That is, you can assume equality of variances if the ratio of the
largest sample standard deviation to the smallest is less than 2. For the chick
weights data, the following code will determine this:
R> chick.sds <- tapply(chickwts$weight,INDEX=chickwts$feed,FUN=sd)
R> chick.sds
casein horsebean linseed meatmeal soybean sunflower
64.43384 38.62584 52.23570 64.90062 54.12907 48.83638
R> max(chick.sds)/min(chick.sds)
[1] 1.680238
Analysis of Variance 437
This informal result indicates that it's reasonable to make the
assumption.
Next, consider the assumption of normality of the raw observations.
This can be difficult to determine in many real-data examples. At the least,
though, it's worthwhile to inspect histograms and QQ plots for signs of nonnormality. You already inspected histograms and QQ plots for all 71 weights
in Section 16.2.2, but for an ANOVA, you need to do this with respect to the
grouping of the observations (that is, not just "overall" for the whole set of
weights regardless of groups).
To achieve this for the chickwts data, you need to first mean-center each
weight by its respective sample mean. You can do this by taking the original vector of weights and subtracting from it the chick.means vector, but
first you must rearrange and replicate the latter elements to correspond to
the elements in the former. This is done by using as.numeric on the factor
vector that represents feed type, giving the numeric value of the levels of
chickwts$feed for each record in the original data frame. When that numeric
vector is passed via the square brackets to chick.means, you get the correct
group mean matched to each observation. As an exercise, you can inspect
all the ingredients that go into creating the following chick.meancen object to
satisfy yourself of what's going on:
R> chick.meancen <- chickwts$weight-chick.means[as.numeric(chickwts$feed)]
In the context of the current analysis, these group-wise, mean-centered
values are also referred to as residuals, a term you'll come across frequently
when you study regression methods in the next few chapters.
You can now assess normality of the observations as a whole using the
residuals. To inspect a normal QQ plot, the relevant functions are qqnorm
and qqline, which you first met in Section 16.2.2. The following two lines
produce the image on the right of Figure 19-1.
R> qqnorm(chick.meancen,main="Normal QQ plot of residuals")
R> qqline(chick.meancen)
Based on this plot (the proximity of the plotted points to the perfect
straight line), it doesn't seem unreasonable to assume normality for these
data, particularly when compared to QQ plots of generated normal data of
the same sample size (an example was given on the left of Figure 16-9 on
page 355).
Investigating the validity of any required assumptions is referred to as
diagnostic checking. If you wanted to perform a more rigorous diagnostic
check for an ANOVA, other visual diagnostics could involve inspecting QQ
plots split by group (you'll do this in an example in Section 19.3) or plotting
the sample standard deviation for each group against the corresponding
sample means. Indeed, there are also general hypothesis tests for normality
(such as the Shapiro-Wilk test or Anderson-Darling test-you'll see the former used in Section 22.3.2), as well as tests for equality of variances (such as
438 Chapter 19
Levene's test), but I'll stick with the basic rule of thumb and visual checks in
this example.
19.1.2 One-Way ANOVA Table Construction
Turning your attention back to the left of Figure 19-1, remember that the
goal is to statistically evaluate the equality of the means marked by ×. This
task will therefore require you to consider not only the variability within each
of the k samples but the variability between the samples; this is why the test is
referred to as an analysis of variance.
The test proceeds by first calculating various metrics associated with
the overall variability and then calculating the within- and between-group
variability. These figures involve sums of squared quantities and associated
degrees of freedom values. All this culminates in a single test statistic and
p-value targeting the aforementioned hypotheses. These ingredients are
typically presented in a table, which is defined as follows.
Let x1, . . . , xN represent all N observations, regardless of group; let
x1(j)
, . . . , xnj (j) denote the specific group observations in group j = 1, . . . , k
such that n1 + . . . + nk = N. Let the "grand mean" of all observations be
defined as x¯ =
1
N
PN
i=1
xi
. The ANOVA table is then constructed, where SS
stands for sum-of-squares, df stands for degrees of freedom, MS stands for
mean square, F refers to the F test statistic, and p refers to the p-value.
df SS MS F p
Overall 1 (1)
Group (or "Factor") k ??? 1 (2) (5) (5)÷(6) p-value
Error (or "Residual") N ??? k (3) (6)
TOTAL N (4)
You calculate the values with these formulas:
(1): Nx¯
2
(2):
Pk
j=1
Pn j
i=1
xi(j)
2
nj
(3): (4)???(2)???(1)
(4):
PN
i=1
x
2
i
(5): (2)÷(k ??? 1)
(6): (3)÷(N ??? k)
There are three input sources that are assumed to make up the observed
data, which, when added together, result in the TOTAL row. Let's think
about these in a little more detail:
Overall row This relates to the scale on which the data as a whole sit.
It doesn't affect the outcome of the hypothesis test (since you're interested only in the relative differences between means) and is sometimes
removed from the table, affecting the TOTAL values accordingly.
Analysis of Variance 439
Group row/Factor row This relates to the data in the individual groups
of interest, thereby accounting for the between-group variability.
Error row/Residual row This accounts for the random deviation from
the estimated means of each group, thereby accounting for the withingroup variability.
TOTAL row This represents the raw data, based on the previous three
ingredients. It is used to find the Error SS by differencing.
The three input sources each have a corresponding degrees of freedom
(df) value in the first column and a sum-of-squares (SS) value attached to
the df in the second column. Between- and within-group variability is averaged by dividing the SS by the df, giving the mean squared (MS) component for these two items. The test statistic, F, is found by dividing the mean
squared group (MSG) effect by the mean squared error (MSE) effect. This
test statistic follows the F-distribution (refer to Section 16.2.5), which itself
requires a pair of degrees of freedom values ordered as df1 (which represents the Group df, k???1) and df2 (which represents the Error df, N???k). Like
the chi-squared distribution, the F-distribution is unidirectional in nature,
and the p-value is obtained as the upper-tail area from the test statistic F.
19.1.3 Building ANOVA Tables with the aov Function
As you might expect, R allows you to easily construct an ANOVA table for
the chick weight test using the built-in aov function as follows:
R> chick.anova <- aov(weight~feed,data=chickwts)
Then, the table is printed to the console screen using summary.
R> summary(chick.anova)
Df Sum Sq Mean Sq F value Pr(>F)
feed 5 231129 46226 15.37 5.94e-10 ***
Residuals 65 195556 3009
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
There are several comments to make here. Note that you employ formula notation weight~feed to specify the measurement variable of interest,
weight, as modeled by the categorical-nominal variable of interest, feed type.
In this case, the variable names weight and feed are not required to be prefaced by chickwts$ since the optional data argument has been passed the data
frame of interest.
Remember from Section 14.3.2 that for the notation in the expression
weight~feed, the "outcome" of interest must always appear on the left of the ~
(this notation will become particularly relevant in Chapters 20 to 22).
440 Chapter 19
To actually view the table, you must apply the summary command to the
object resulting from the call to aov. R omits the first and last rows (Overall
and TOTAL) since these are not directly involved in calculating the p-value.
Other than that, it's easy to identify that the feed row refers to the Group row
and the Residuals row refers to the Error row.
NOTE By default, R annotates model-based summary output like this with significance stars.
These show intervals of significance, and the number of stars increases as the p-value
decreases beyond a cutoff mark of 0.1. This can be useful when you're examining more
complicated analyses where multiple p-values are summarized, though not everyone
likes this feature. If you want, you can turn off this feature in a given R session by
entering options(show.signif.stars=FALSE) at the prompt. Alternatively, you can
turn off the feature directly in the call to summary by setting the additional argument
signif.stars=FALSE. In this book, I'll leave them be.
From the contents of the ANOVA for this example, you can quickly confirm the calculations. Note that the MSE, 3009, was defined as the Error SS
divided by the Error df. Indeed, in R, the same result is achieved manually
(the table output has been rounded to the nearest integer).
R> 195556/65
[1] 3008.554
You can confirm all the other results in the table output using the relevant equations from earlier.
Interpreting a hypothesis test based on ANOVA follows the same rules
as any other test. With the understanding of a p-value as "the probability
that you observe the sample statistics at hand or something more extreme,
if H0 is true," a small p-value indicates evidence against the null hypothesis. In the current example, a tiny p-value provides strong evidence against
the null that the mean chick weights are the same for the different diets.
In other words, you reject H0 in favor of HA; the latter states that there is a
difference.
In a similar fashion as in the chi-squared tests, rejection of the null in
one-way ANOVA doesn't tell you exactly where a difference lies, merely that
there's evidence one exists. Further scrutiny of the data in the individual
groups is necessary to identify the offending means. At the simplest level,
you could turn back to pairwise two-sample t-tests, in which case you could
also use the MSE from the ANOVA table as an estimate of the pooled variance. The substitution is valid if the assumption of equal variance holds, and
such a step is beneficial because the corresponding t-based sampling distribution will utilize the Error df (this is naturally higher than would otherwise
be the case if the df was based on just the sample sizes of the two groups of
specific interest).
Analysis of Variance 441
Exercise 19.1
Consider the following data:
Site I Site II Site III Site IV
93 85 100 96
120 45 75 58
65 80 65 95
105 28 40 90
115 75 73 65
82 70 65 80
99 65 50 85
87 55 30 95
100 50 45 82
90 40 50
78 45
95 55
93
88
110
These figures provide the depths (in centimeters) at which
important archaeological finds were made at four sites in New Mexico (see Woosley and Mcintyre, 1996). Store these data in your R
workspace, with one vector containing depth and the other vector
containing the site of each observation.
a. Produce side-by-side boxplots of the depths split by group, and
use additional points to mark the locations of the sample means.
b. Assuming independence, execute diagnostic checks for normality and equality of variances.
c. Perform and conclude a one-way ANOVA test for evidence of a
difference between the means.
In Section 14.4, you looked at the data set providing measurements
on petal and sepal sizes for three species of iris flowers. This is available in R as iris.
d. Based on diagnostic checks for normality and equality of variances, decide which of the four outcome measurements (sepal
length/width and petal length/width) would be suitable for
ANOVA (using the species as the group variable).
e. Carry out one-way ANOVA for any suitable measurement
variables.
442 Chapter 19
19.2 Two-Way ANOVA
In many studies, the numeric outcome variable you're interested in will
be categorized by more than just one grouping variable. In these cases,
you would use the multiple-factor ANOVA rather than the one-way ANOVA.
This technique is directly referred to by the number of grouping variables
used, with two- and three-way ANOVA being the next and most common
extensions.
Increasing the number of grouping variables complicates matters
somewhat-performing just a one-way ANOVA for each variable separately is inadequate. In dealing with more than one categorical grouping
factor, you must consider the main effects of each factor on the numeric outcome, while simultaneously accounting for the presence of the other grouping factor(s). That's not all, though. It's just as important to additionally
investigate the idea of an interactive effect; if an interactive effect exists, then
it suggests that the impact one of the grouping variables has on the outcome
of interest, specified by its main effect, varies according to the levels of the
other grouping variable(s).
19.2.1 A Suite of Hypotheses
For this explanation, denote your numeric outcome variable with O and
your two grouping variables as G1 and G2. In two-way ANOVA, the hypotheses should be set along the following lines:
H0 : G1 has no main (marginal) effect on the mean of O.
G2 has no main (marginal) effect on the mean of O.
There is no interactive effect of G1 with G2 on the mean of O.
HA : Separately, each statement in H0 is incorrect.
You can see from these general hypotheses that you now have to obtain
a p-value for each of the three components.
For the example, let's use the built-in warpbreaks data frame (Tippett,
1950), which provides the number of "warp break" imperfections (column
breaks) observed in 54 pieces of yarn of equal length. Each piece of yarn is
classified according to two categorical variables: wool (the type of yarn, with
levels A and B) and tension (the level of tension applied to that piece-L, M, or
H for low, medium, or high). Using tapply, you can inspect the mean number
of warp breaks for each classification.
R> tapply(warpbreaks$breaks,INDEX=list(warpbreaks$wool,warpbreaks$tension),
FUN=mean)
L M H
A 44.55556 24.00000 24.55556
B 28.22222 28.77778 18.77778
Analysis of Variance 443
You can supply more than one grouping variable to the INDEX argument as separate members of a list (any factor vectors given to this argument
should be the same length as the first argument that specifies the data of
interest). The results are returned as a matrix for two grouping variables, a
3D array for three grouping variables, and so on.
For some analyses, however, you might need the same information provided earlier in a different format. The aggregate function is similar to tapply,
but it returns a data frame, the results in stacked format according to the
specified grouping variables (as opposed to an array as returned by tapply).
It's called in much the same way. The first argument takes the data vector of
interest. The second argument, by, should be a list of the desired grouping
variables, and in FUN, you specify the function to operate on each subset.
R> wb.means <- aggregate(warpbreaks$breaks,
by=list(warpbreaks$wool,warpbreaks$tension),FUN=mean)
R> wb.means
Group.1 Group.2 x
1 A L 44.55556
2 B L 28.22222
3 A M 24.00000
4 B M 28.77778
5 A H 24.55556
6 B H 18.77778
Here I've stored the result of the call to aggregate as the object wb.means
for later use.
19.2.2 Main Effects and Interactions
I mentioned earlier that you could perform just a one-way ANOVA on each
grouping variable separately, but this, in general, isn't a good idea. I'll
demonstrate this now with the warpbreaks data (a quick inspection of the
relevant diagnostics shows no obvious cause for concern):
R> summary(aov(breaks~wool,data=warpbreaks))
Df Sum Sq Mean Sq F value Pr(>F)
wool 1 451 450.7 2.668 0.108
Residuals 52 8782 168.9
R> summary(aov(breaks~tension,data=warpbreaks))
Df Sum Sq Mean Sq F value Pr(>F)
tension 2 2034 1017.1 7.206 0.00175 **
Residuals 51 7199 141.1
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
This output tells you that if you ignore tension, there is no evidence to
suggest that there is any difference in the mean number of imperfections
444 Chapter 19
based on the type of wool alone (p-value 0.108). If you ignore wool, however, there is evidence to suggest a difference in warp breaks according to
tension only.
The problem here is that by ignoring one of the variables, you lose the
ability to detect differences (or, more generally, statistical relationships) that
may occur at a finer level. For example, though the wool type alone seems to
have no remarkable impact on the mean number of warp breaks, you cannot tell whether this would be the case if you just looked at wool types at one
particular level of tension.
Instead, you investigate this kind of question using two-way ANOVA. The
following executes a two-way ANOVA for the warp breaks data based only on
the main effects of the two grouping variables:
R> summary(aov(breaks~wool+tension,data=warpbreaks))
Df Sum Sq Mean Sq F value Pr(>F)
wool 1 451 450.7 3.339 0.07361 .
tension 2 2034 1017.1 7.537 0.00138 **
Residuals 50 6748 135.0
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Take a look at the formula. Specifying wool+tension to the right of the
outcome variable and the ~ allows you to take both grouping variables into
account at the same time. The results reveal a small drop in the size of
the p-values now attached to each grouping variable; indeed, the p-value
for wool is around 0.073, approaching the conventional cutoff significance
level of ?? = 0.05. To interpret the results, you hold one grouping variable
constant-if you focus on just one type of wool, there is still statistically significant evidence to suggest a difference in the mean number of warp breaks
between the different tension levels. If you focus on just one level of tension,
the evidence of a difference considering the two wool types has increased a
little but is still not statistically significant (assuming the aforementioned
?? = 0.05).
There's still a limitation with considering only main effects. While the
previous analysis shows that there's variation in the outcome between the
different levels of the two categorical variables, it doesn't address the possibility that a difference in the mean number of warp breaks might vary
further according to precisely which level of either tension or wool is being
used when holding the other variable constant. This relatively subtle yet
important consideration is known as an interaction. Specifically, if there is
an interactive effect present between tension and wool with respect to warp
breaks, then this would imply that the magnitude and/or direction of the
difference in the mean number of warp breaks is not the same at different
levels of the two grouping factors.
To account for interactions, you make a slight adjustment to the two-way
ANOVA model code.
Analysis of Variance 445
R> summary(aov(breaks~wool+tension+wool:tension,data=warpbreaks))
Df Sum Sq Mean Sq F value Pr(>F)
wool 1 451 450.7 3.765 0.058213 .
tension 2 2034 1017.1 8.498 0.000693 ***
wool:tension 2 1003 501.4 4.189 0.021044 *
Residuals 48 5745 119.7
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
You can explicitly specify the interaction as the main effects model
formula plus the notation wool:tension, where the two grouping variables
are separated by a :. (Note, in this setting, the : operator has nothing to
do with the shortcut for creating an integer sequence first discussed in
Section 2.3.2.)
You can see from the ANOVA table output that, statistically, there is evidence of an interactive effect; that is, the very nature of the difference in
the means is dependent upon the factor levels themselves, even though that
evidence is relatively weak. Of course, the p-value of around 0.021 tells you
only that, overall, there might be an interaction but not the precise features
of the interaction.
To help with this, you can interpret such a two-way interaction effect in
more detail with an interaction plot, provided in R with interaction.plot.
R> interaction.plot(x.factor=wb.means[,2],trace.factor=wb.means[,1],
response=wb.means$x,trace.label="wool",
xlab="tension",ylab="mean warp breaks")
When interaction.plot is called, the outcome means should be supplied
to the argument response, and the vectors providing the corresponding levels
of each of the two factors should be supplied to the arguments x.factor (for
the variable on the horizontal axis that refers to moving between levels from
the left to the right) and trace.factor (each level of which will produce a different line, referenced in an automatically produced legend; the title of this
legend is passed to trace.label). It doesn't matter which grouping variable is
which; the appearance of the plot will change accordingly, but your interpretations will (should!) be the same. The result is shown in Figure 19-2.
The two-way interaction plot displays the outcome variable on the vertical axis and splits the recorded means by the levels of the two grouping
variables. This allows you to inspect the potential effect that varying the
levels of the grouping variables has on the outcome. In general, when the
lines (or segments thereof) are not parallel, it suggests an interaction could
be present. Vertical separations between the plotted locations indicate the
individual main effects of the grouping variables.
It turns out that the columns returned by a call to aggregate are actually
perfectly suited to interaction.plot. As usual, you can specify the common
graphical parameters, like those you initially encountered in Section 7.2, to
control specific features of the plot and axis annotation. For Figure 19-2,
446 Chapter 19
you've specified that x.factor should be the second column of the wb.means
matrix, meaning that the tension levels vary horizontally. The trace.factor
here is the type of wool, so there are only two distinct lines corresponding to
the two levels A and B. The response is that third column of wb.means, extracted
using $x (take a look at the wb.means object; you'll see the column containing
the results of interest is labeled x by default after a call to aggregate).
Figure 19-2: An interaction plot for the full two-way ANOVA
model of the warpbreaks data set
Considering the actual appearance of the plot in Figure 19-2, it does
indeed appear that the mean number of warp breaks for wool type A is
higher if tension is low, but the nature of the difference changes if you move
to a medium tension, where B has a higher point estimate than A. Moving to
a high tension, type A again has a higher estimate of the mean number of
breaks than B, though here the difference between A and B is nowhere near
as big as it is at a low tension. (Note, however, that the interaction plot does
not display any kind of standard error measurements, so you must remember that all point estimates of the means are subject to variability.)
Interactions are certainly not a concept unique to multiple-factor
ANOVA; they form an important consideration in many different types of
statistical models. For the moment, it's good just to gain a basic appreciation
of interactions.
19.3 Kruskal-Wallis Test
When comparing multiple means, there may be situations when you're
unwilling to assume normality or have even found the assumption of normality invalid in diagnostic checks. In this case, you can use the KruskalWallis test, an alternative to the one-way ANOVA that relaxes the dependence
Analysis of Variance 447
on the necessity for normality. This method tests for "equality of distributions" of the measurements in each level of the grouping factor. If you make
the usual assumption of equal variances across these groups, you can therefore think of this test as one that compares multiple medians rather than
means.
The hypotheses governing the test alter accordingly.
H0 : Group medians are all equal.
HA : Group medians are not all equal
(alternatively, at least one group median differs).
The Kruskal-Wallis test is a nonparametric approach since it does not
rely on quantiles of a standardized parametric distribution (in other
words, the normal distribution) or any of its functions. In the same way
that the ANOVA is a generalization of the two-sample t-test, the KruskalWallis ANOVA is a generalization of the Mann-Whitney test for two medians.
It's also referred to as the Kruskal-Wallis rank sum test, and you use the chisquared distribution to calculate the p-value.
Turn your attention to the data frame survey, located in the MASS package. These data record particular characteristics of 237 first-year undergraduate statistics students collected from a class at the University of Adelaide,
South Australia. Load the required package first with a call to library("MASS")
and then enter ?survey at the prompt. You can read the help file to understand which variables are present in the data frame.
Suppose you're interested to see whether the age of the students, Age,
tends to differ with respect to four smoking categories reported in Smoke.
An inspection of the relevant side-by-side boxplots and a normal QQ plot
of the residuals (mean-centered observations with respect to each group)
suggests a straightforward one-way ANOVA isn't necessarily a good idea. The
following code (which mimics the steps you saw in Section 19.1.1) produces
the two images in Figure 19-3, which show normality is questionable:
R> boxplot(Age~Smoke,data=survey)
R> age.means <- tapply(survey$Age,survey$Smoke,mean)
R> age.meancen <- survey$Age-age.means[as.numeric(survey$Smoke)]
R> qqnorm(age.meancen,main="Normal QQ plot of residuals")
R> qqline(age.meancen)
With this possible violation of normality, you could therefore apply the
Kruskal-Wallis test instead of the parametric ANOVA. A quick check for
equality of variances further supports this, with the ratio of the largest to
the smallest group standard deviations clearly being less than 2.
R> tapply(survey$Age,survey$Smoke,sd)
Heavy Never Occas Regul
6.332628 6.675257 5.861992 5.408822
448 Chapter 19
Figure 19-3: Side-by-side boxplots (left) and a normal QQ plot of the residuals (right) for
the student age observations split by smoking status
In R, a Kruskal-Wallis test is performed using kruskal.test.
R> kruskal.test(Age~Smoke,data=survey)
Kruskal-Wallis rank sum test
data: Age by Smoke
Kruskal-Wallis chi-squared = 3.9262, df = 3, p-value = 0.2695
The syntax for this test is the same as for aov. As you might suspect from
Figure 19-3, the large p-value suggests there's no evidence against the null
hypothesis that states that the medians are all equal. In other words, there
doesn't seem to be an overall age difference between the students in the
four smoking categories.
Exercise 19.2
Bring up the quakes data frame again, which describes the locations, magnitudes, depths, and number of observation stations that
detected 1,000 seismic events off the coast of Fiji.
a. Use cut (see Section 4.3.3) to create a new factor vector defining the depths of each event according to the following three
categories: (0,200], (200,400], and (400,680].
b. Decide whether a one-way ANOVA or a Kruskal-Wallis test is
more appropriate to use to compare the distributions of the
number of detecting stations, split according to the three categories in (a).
Analysis of Variance 449
c. Perform your choice of test in (b) (assume a ?? = 0.01 level of significance) and conclude.
Load the MASS package with a call to library("MASS") if you haven't
already done so in the current R session. This package includes
the ready-to-use Cars93 data frame, which contains detailed data on
93 cars for sale in the United States in 1993 (Lock, 1993; Venables
and Ripley, 2002).
d. Use aggregate to compute the mean length of the 93 cars, split by
two categorical variables: AirBags (type of airbags available-levels
are Driver & Passenger, Driver only, and None), and Man.trans.avail
(whether the car comes in a manual transmission-levels are Yes
and No).
e. Produce an interaction plot using the results in (d). Does there
appear to be an interactive effect of AirBags with Man.trans.avail
on the mean length of these cars (if you consider only these
variables)?
f. Fit a full two-way ANOVA model for the mean lengths according
to the two grouping variables (assume satisfaction of all relevant
assumptions). Is the interactive effect statistically significant? Is
there evidence of any main effects?
Important Code in This Chapter
Function/operator Brief description First occurrence
aov Produce ANOVA table Section 19.1.3, p. 440
aggregate Stacked statistics by factor Section 19.2.1, p. 444
interaction.plot Two-factor interaction plot Section 19.2.2, p. 446
kruskal.test Kruskal-Wallis test Section 19.3, p. 449
450 Chapter 19
20
S IMPLE L INEAR REGRESS ION
Though straightforward comparative tests
of individual statistics are useful in their
own right, you'll often want to learn more
from your data. In this chapter, you'll look at
linear regression models: a suite of methods used to evaluate precisely how variables relate to each other.
Simple linear regression models describe the effect that a particular variable, called the explanatory variable, might have on the value of a continuous
outcome variable, called the response variable. The explanatory variable may
be continuous, discrete, or categorical, but to introduce the key concepts,
I'll concentrate on continuous explanatory variables for the first several sections in this chapter. Then, I'll cover how the representation of the model
changes if the explanatory variable is categorical.
20.1 An Example of a Linear Relationship
As an example to start with, let's continue with the data used in Section 19.3
and look at the student survey data (the survey data frame in the package
MASS) a little more closely. If you haven't already done so, with the required
package loaded (call library("MASS")), you can read the help file ?survey for
details on the variables present.
Plot the student heights on the y-axis and their handspans (of their writing hand) on the x-axis.
R> plot(survey$Height~survey$Wr.Hnd,xlab="Writing handspan (cm)",
ylab="Height (cm)")
Figure 20-1 shows the result.
Figure 20-1: A scatterplot of height against writing handspan
for a sample of first-year statistics students
Note that the call to plot uses formula notation (also referred to as
symbolic notation) to specify "height on handspan." You can produce the
same scatterplot by using the coordinate vector form of (x, y), that is,
plot(survey$Wr.Hnd,survey$Height,...), but I'm using the symbolic notation
here because it nicely reflects how you'll fit the linear model in a moment.
As you might expect, there's a positive association between a student's
handspan and their height. That relationship appears to be linear in nature.
To assess the strength of the linear relationship (refer to Section 13.2.5), you
can find the estimated correlation coefficient.
R> cor(survey$Wr.Hnd,survey$Height,use="complete.obs")
[1] 0.6009909
Though there are 237 records in the data frame, the plot doesn't
actually show 237 points. This is because there are missing observations
(coded NA; see Section 6.1.3). By default, R removes any "incomplete" pairs
when producing a plot like this. To find out how many offending observations have been deleted, you can use the short-form logical operator |
452 Chapter 20
(Section 4.1.3) in conjunction with is.na (Section 6.1.3) and which (Section 4.1.5). You then use length to discover there are 29 missing observation
pairs.
R> incomplete.obs <- which(is.na(survey$Height)|is.na(survey$Wr.Hnd))
R> length(incomplete.obs)
[1] 29
NOTE Because there are NAs in the vectors supplied to the correlation coefficient function cor,
you must also specify the optional argument use="complete.obs". This means that the
calculated statistic takes into account only those observation pairs in the Wr.Hnd and
Height vectors where neither element is NA. You can think of this argument as doing
much the same thing as na.rm=TRUE in univariate summary statistic functions such as
mean and sd.
20.2 General Concepts
The purpose of a linear regression model is to come up with a function
that estimates the mean of one variable given a particular value of another
variable. These variables are known as the response variable (the "outcome"
variable whose mean you are attempting to find) and the explanatory variable
(the "predictor" variable whose value you already have).
In terms of the student survey example, you might ask something like
"What's the expected height of a student if their handspan is 14.5 cm?"
Here the response variable is the height, and the explanatory variable is the
handspan.
20.2.1 Definition of the Model
Assume you're looking to determine the value of response variable Y given
the value of an explanatory variable X. The simple linear regression model states
that the value of a response is expressed as the following equation:
Y |X = ??0 + ??1X + o (20.1)
On the left side of Equation (20.1), the notation Y |X reads as "the value
of Y conditional upon the value of X."
Residual Assumptions
The validity of the conclusions you can draw based on the model in (20.1) is
critically dependent on the assumptions made about o, which are defined as
follows:
. The value of o is assumed to be normally distributed such that
o ??? N(0,??).
. That o is centered (that is, has a mean of) zero.
. The variance of o, ??
2
, is constant.
Simple Linear Regression 453
The o term represents random error. In other words, you assume that
any raw value of the response is owed to a linear change in a given value
of X, plus or minus some random, residual variation or normally distributed
noise.
Parameters
The value denoted by ??0 is called the intercept, and that of ??1 is called the
slope. Together, they are also referred to as the regression coefficients and are
interpreted as follows:
. The intercept, ??0, is interpreted as the expected value of the response
variable when the predictor is zero.
. Generally, the slope, ??1, is the focus of interest. This is interpreted as
the change in the mean response for each one-unit increase in the predictor. When the slope is positive, the regression line increases from left
to right (the mean response is higher when the predictor is higher);
when the slope is negative, the line decreases from left to right (the
mean response is lower when the predictor is higher). When the slope
is zero, this implies that the predictor has no effect on the value of the
response. The more extreme the value of ??1 (that is, away from zero),
the steeper the increasing or decreasing line becomes.
20.2.2 Estimating the Intercept and Slope Parameters
The goal is to use your data to estimate the regression parameters, yielding
the estimates ??^
0 and ??^
1; this is referred to as fitting the linear model. In this
case, the data comprise n pairs of observations for each individual. The fitted model of interest concerns the mean response value, denoted y^, for a
specific value of the predictor, x, and is written as follows:
y^ = ??^
0 + ??^
1 x (20.2)
Sometimes, alternative notation such as E[Y] or E[Y |X = x] is used on
the left side of (20.2) to emphasize the fact that the model gives the mean
(that is, the expected value) of the response. For compactness, many simply
use something like y^, as shown here.
Let your n observed data pairs be denoted xi and yi for the predictor
and response variables, respectively; i = 1, . . . , n. Then, the parameter estimates for the simple linear regression function are
??^
1 = ??x y
sy
sx
and ??^
0 = y¯ ??? ??^
1 x¯ (20.3)
where
. x¯ and y¯ are the sample means of the xis and yis.
. sx and sy are the sample standard deviations of the xis and yis.
. ??x y is the estimate of correlation between X and Y based on the data
(see Section 13.2.5).
454 Chapter 20
Estimating the model parameters in this way is referred to as least-squares
regression; the reason for this will become clear in a moment.
20.2.3 Fitting Linear Models with lm
In R, the command lm performs the estimation for you. For example, the following line creates a fitted linear model object of the mean student height
by handspan and stores it in your global environment as survfit:
R> survfit <- lm(Height~Wr.Hnd,data=survey)
The first argument is the now-familiar response ~ predictor formula,
which specifies the desired model. You don't have to use the survey$ prefix
to extract the vectors from the data frame because you specifically instruct
lm to look in the object supplied to the data argument.
The fitted linear model object itself, survfit, has a special class in R-
one of "lm". An object of class "lm" can essentially be thought of as a list containing several components that describe the model. You'll look at these in a
moment.
If you simply enter the name of the "lm" object at the prompt, it will
provide the most basic output: a repeat of your call and the estimates of the
intercept (??^
0) and slope (??^
1).
R> survfit
Call:
lm(formula = Height ~ Wr.Hnd, data = survey)
Coefficients:
(Intercept) Wr.Hnd
113.954 3.117
This reveals that the linear model for this example is estimated as
follows:
y^ = 113.954 + 3.117x (20.4)
If you evaluate the mathematical function for y^-Equation (20.2)-at
a range of different values for x, you end up with a straight line when you
plot the results. Considering the definition of intercept given earlier as the
expected value of the response variable when the predictor is zero, in the
current example, this would imply that the mean height of a student with
a handspan of 0 cm is 113.954 cm (an arguably less-than-useful statement
since a value of zero for the explanatory variable doesn't make sense; you'll
consider these and related issues in Section 20.4). The slope, the change
in the mean response for each one-unit increase in the predictor, is 3.117.
This states that, on average, for every 1 cm increase in handspan, a student's
height is estimated to increase by 3.117 cm.
Simple Linear Regression 455
With all this in mind, once more run the line to plot the raw data as
given in Section 20.1 and shown in Figure 20-1, but now add the fitted
regression line using abline. So far, you've only used the abline command
to add perfectly horizontal and vertical lines to an existing plot, but when
passed an object of class "lm" that represents a simple linear model, like
survfit, the fitted regression line will be added instead.
R> abline(survfit,lwd=2)
This adds the slightly thickened diagonally increasing line shown in
Figure 20-2.
Figure 20-2: The simple linear regression line (solid, bold)
fitted to the observed data. Two dashed vertical line
segments provide examples of a positive (leftmost) and
negative (rightmost) residual.
20.2.4 Illustrating Residuals
When the parameters are estimated as shown here, using (20.3), the fitted
line is referred to as an implementation of least-squares regression because it's
the line that minimizes the average squared difference between the observed
data and itself. This concept is easier to understand by drawing the distances
between the observations and the fitted line, formally called residuals, for a
couple of individual observations within Figure 20-2.
First, let's extract two specific records from the Wr.Hnd and Height data
vectors and call the resulting vectors obsA and obsB.
R> obsA <- c(survey$Wr.Hnd[197],survey$Height[197])
R> obsA
[1] 15.00 170.18
456 Chapter 20
R> obsB <- c(survey$Wr.Hnd[154],survey$Height[154])
R> obsB
[1] 21.50 172.72
Next, briefly inspect the names of the members of the survfit object.
R> names(survfit)
[1] "coefficients" "residuals" "effects" "rank"
[5] "fitted.values" "assign" "qr" "df.residual"
[9] "na.action" "xlevels" "call" "terms"
[13] "model"
These members are the components that automatically make up a fitted model object of class "lm", mentioned briefly earlier. Note that there's
a component called "coefficients". This contains a numeric vector of the
estimates of the intercept and slope.
You can extract this component (and indeed any of the other ones listed
here) in the same way you would perform a member reference on a named
list: by entering survfit$coefficients at the prompt. Where possible, though,
it's technically preferable for programming purposes to extract such components using a "direct-access" function. For the coefficients component of an
"lm" object, the function you use is coef.
R> mycoefs <- coef(survfit)
R> mycoefs
(Intercept) Wr.Hnd
113.953623 3.116617
R> beta0.hat <- mycoefs[1]
R> beta1.hat <- mycoefs[2]
Here, the regression coefficients are extracted from the object and then
separately assigned to the objects beta0.hat and beta1.hat. Other common
direct-access functions include resid and fitted; these two pertain to the
"residuals" and "fitted.values" components, respectively.
Finally, I use segments to draw the vertical dashed lines present in
Figure 20-2.
R> segments(x0=c(obsA[1],obsB[1]),y0=beta0.hat+beta1.hat*c(obsA[1],obsB[1]),
x1=c(obsA[1],obsB[1]),y1=c(obsA[2],obsB[2]),lty=2)
Note that the dashed lines meet the fitted line at the vertical axis locations passed to y0, which, with the use of the regression coefficients beta0.hat
and beta1.hat, reflects Equation (20.4).
Now, imagine a collection of alternative regression lines drawn through
the data (achieved by altering the value of the intercept and slope). Then,
for each of the alternative regression lines, imagine you calculate the residuals (vertical distances) between the response value of every observation and
Simple Linear Regression 457
the fitted value of that line. The simple linear regression line estimated as
per (20.3) is the line that lies "closest to all observations." By this, it is meant
that the fitted regression model is represented by the estimated line that
passes through the coordinate provided by the variable means (x¯, y¯), and
it's the line that yields the smallest overall measure of the squared residual
distances. For this reason, another name for a least-squares-estimated regression equation like this is the line of best fit.
20.3 Statistical Inference
The estimation of a regression equation is relatively straightforward, but this
is merely the beginning. You should now think about what can be inferred
from your result. In simple linear regression, there's a natural question that
should always be asked: Is there statistical evidence to support the presence
of a relationship between the predictor and the response? To put it another
way, is there evidence that a change in the explanatory variable affects the
mean outcome? You investigate this following the same ideas that were introduced in Chapter 17 when you began thinking about the variability present
in estimated statistics and then continued to infer from your results using
confidence intervals and, in Chapter 18, hypothesis testing.
20.3.1 Summarizing the Fitted Model
This kind of model-based inference is automatically carried out by R when lm
objects are processed. Using the summary function on an object created by lm
provides you with output far more detailed than simply printing the object
to the console. For the moment, you'll focus on just two aspects of the information presented in summary: the significance tests associated with the regression coefficients and the interpretation of the so-called coefficient of determination (labeled R-squared in the output), which I'll explain shortly.
Use summary on the current model object survfit, and you'll see the
following:
R> summary(survfit)
Call:
lm(formula = Height ~ Wr.Hnd, data = survey)
Residuals:
Min 1Q Median 3Q Max
-19.7276 -5.0706 -0.8269 4.9473 25.8704
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 113.9536 5.4416 20.94 <2e-16 ***
Wr.Hnd 3.1166 0.2888 10.79 <2e-16 ***
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
458 Chapter 20
Residual standard error: 7.909 on 206 degrees of freedom
(29 observations deleted due to missingness)
Multiple R-squared: 0.3612, Adjusted R-squared: 0.3581
F-statistic: 116.5 on 1 and 206 DF, p-value: < 2.2e-16
20.3.2 Regression Coefficient Significance Tests
Let's begin by focusing on the way the estimated regression coefficients are
reported. The first column of the coefficients table contains the point estimates of the intercept and slope (the intercept is labeled as such, and the
slope is labeled after the name of the predictor variable in the data frame);
the table also includes estimates of the standard errors of these statistics.
It can be shown that simple linear regression coefficients, when estimated
using least-squares, follow a t-distribution with n ??? 2 degrees of freedom
(when given the number of observations, n, used in the model fit). The
standardized t value and a p-value are reported for each parameter. These
represent the results of a two-tailed hypothesis test formally defined as
H0 : ??j = 0
HA : ??j , 0
where j = 0 for the intercept and j = 1 for the slope, using the notation in
Equation (20.1).
Focus on the row of results for the predictor. With a null value of zero,
truth of H0 implies that the predictor has no effect on the response. The
claim here is interested in whether there is any effect of the covariate, not
the direction of this effect, so HA is two-sided (via ,). As with any hypothesis
test, the smaller the p-value, the stronger the evidence against H0. With a
small p-value (< 2 × 10???16) attached to this particular test statistic (which
you can confirm using the formula in Chapter 18: T = (3.116 ??? 0)/0.2888 =
10.79), you'd therefore conclude there is strong evidence against the claim
that the predictor has no effect on the mean level of the response.
The same test is carried out for the intercept, but the test for the
slope parameter ??1 is typically more interesting (since rejection of the null
hypothesis for ??0 simply indicates evidence that the regression line does
not strike the vertical axis at zero), especially when the observed data don't
include x = 0, as is the case here.
From this, you can conclude that the fitted model suggests there is evidence that an increase in handspan is associated with an increase in height
among the population being studied. For each additional centimeter of
handspan, the average increase in height is approximately 3.12 cm.
You could also produce confidence intervals for your estimates using
Equation (17.2) on page 378 and knowledge of the sampling distributions
of the regression parameters; however, yet again, R provides a convenient
function for an object of class "lm" to do this for you.
Simple Linear Regression 459
R> confint(survfit,level=0.95)
2.5 % 97.5 %
(Intercept) 103.225178 124.682069
Wr.Hnd 2.547273 3.685961
To the confint function you pass your model object as the first argument
and your desired level of confidence as level. This indicates that you should
be 95 percent confident the true value of ??1 lies somewhere between 2.55
and 3.69 (to 2 d.p.). As usual, the exclusion of the null value of zero reflects
the statistically significant result from earlier.
20.3.3 Coefficient of Determination
The output of summary also provides you with the values of Multiple R-squared
and Adjusted R-squared, which are particularly interesting. Both of these are
referred to as the coefficient of determination; they describe the proportion of
the variation in the response that can be attributed to the predictor.
For simple linear regression, the first (unadjusted) measure is simply
obtained as the square of the estimated correlation coefficient (refer to
Section 13.2.5). For the student height example, first store the estimated
correlation between Wr.Hnd and Height as rho.xy, and then square it:
R> rho.xy <- cor(survey$Wr.Hnd,survey$Height,use="complete.obs")
R> rho.xy^2
[1] 0.3611901
You get the same result as the Multiple R-squared value (usually written
mathematically as R
2
). This tells you that about 36.1 percent of the variation
in the student heights can be attributed to handspan.
The adjusted measure is an alternative estimate that takes into account
the number of parameters that require estimation. The adjusted measure is
generally important only if you're using the coefficient of determination to
assess the overall "quality" of the fitted model in terms of a balance between
goodness of fit and complexity. I'll cover this in Chapter 22, so I won't go
into any more detail just yet.
20.3.4 Other summary Output
The summary of the model object provides you with even more useful information. The "residual standard error" is the estimated standard error of the o
term (in other words, the square root of the estimated variance of o, namely,
??
2
); below that it also reports any missing values. (The 29 observation pairs
"deleted due to missingness" here matches the number of incomplete observations determined in Section 20.1.)
The output also provides a five-number summary (Section 13.2.3) for
the residual distances-I'll cover this further in Section 22.3. As the final
460 Chapter 20
result, you're provided with a certain hypothesis test performed using the
F -distribution. This is a global test of the impact of your predictor(s) on
the response; this will be explored alongside multiple linear regression in
Section 21.3.5.
You can access all the output provided by summary directly, as individual
R objects, rather than having to read them off the screen from the entire
printed summary. Just as names(survfit) provides you with an indication
of the contents of the stand-alone survfit object, the following code gives
you the names of all the components accessible after summary is used to process survfit.
R> names(summary(survfit))
[1] "call" "terms" "residuals" "coefficients"
[5] "aliased" "sigma" "df" "r.squared"
[9] "adj.r.squared" "fstatistic" "cov.unscaled" "na.action"
It's fairly easy to match most of the components with the printed summary
output, and they can be extracted using the dollar operator as usual. The
residual standard error, for example, can be retrieved directly with this:
R> summary(survfit)$sigma
[1] 7.90878
There are further details on this in the ?summary.lm help file.
20.4 Prediction
To wrap up these preliminary details of linear regression, you'll now look at
using your fitted model for predictive purposes. The ability to fit a statistical
model means that you not only can understand and quantify the nature of
relationships in your data (like the estimated 3.1166 cm increase in mean
height per 1 cm increase in handspan for the student example) but can also
predict values of the outcome of interest, even where you haven't actually
observed the values of any explanatory variables in the original data set. As
with any statistic, though, there is always a need to accompany any point estimates or predictions with a measure of spread.
20.4.1 Confidence Interval or Prediction Interval?
With a fitted simple linear model you're able to calculate a point estimate of
the mean response value, conditional upon the value of an explanatory variable. To do this, you simply plug in (to the fitted model equation) the value
of x you're interested in. A statistic like this is always subject to variation, so
just as with sample statistics explored in earlier chapters, you use a confidence
interval for the mean response (CI) to gauge this uncertainty.
Simple Linear Regression 461
Assume a simple linear regression line has been fitted to n observations
such that y^ = ??^
0+ ??^
1 x. A 100(1?????) percent confidence interval for the mean
response given a value of x is calculated with
y^ ± t(1?????/2,n???2)so
s
1
n
+
(x ??? x¯)
2
(n ??? 1)s
2
x
(20.5)
where you obtain the lower limit by subtraction, the upper limit by addition.
Here, y^ is the fitted value (from the regression line) at x; t(1?????/2,n???2)
is
the appropriate critical value from a t-distribution with n ??? 2 degrees of freedom (in other words, resulting in an upper-tail area of exactly ??/2); so is the
estimated residual standard error; and x¯ and s
2
x
represent the sample mean
and the variance of the observations of the predictor, respectively.
A prediction interval (PI) for an observed response is different from the
confidence interval in terms of context. Where CIs are used to describe the
variability of the mean response, a PI is used to provide the possible range of
values that an individual realization of the response variable might take, given
x. This distinction is subtle but important: the CI corresponds to a mean,
and the PI corresponds to an individual observation.
Let's remain with the previous notation. It can be shown that 100(1 ??? ??)
percent prediction interval for an individual response given a value of x is
calculated with the following:
y^ ± t(1?????/2,n???2)so
s
1 +
1
n
+
(x ??? x¯)
2
(n ??? 1)s
2
x
(20.6)
It turns out that the only difference from (20.5) is the 1+ that appears in
the square root. As such, a PI at x is wider than a CI at x.
20.4.2 Interpreting Intervals
Continuing with our example, let's say you want to determine the mean
height for students with a handspan of 14.5 cm and for students with a
handspan of 24 cm. The point estimates themselves are easy-just plug the
desired x values into the regression equation (20.4).
R> as.numeric(beta0.hat+beta1.hat*14.5)
[1] 159.1446
R> as.numeric(beta0.hat+beta1.hat*24)
[1] 188.7524
According to the model, you can expect mean heights to be around
159.14 and 188.75 cm for handspans of 14.5 and 24 cm, respectively. The
as.numeric coercion function (first encountered in Section 6.2.4) is used
simply to strip the result of the annotative names that are otherwise present
from the beta0.hat and beta1.hat objects.
462 Chapter 20
Confidence Intervals for Mean Heights
To find confidence intervals for these estimates, you could calculate them
manually using (20.5), but of course R has a built-in predict command to do
it for you. To use predict, you first need to store your x values in a particular
way: as a column in a new data frame. The name of the column must match
the predictor used in the original call to create the fitted model object. In
this example, I'll create a new data frame, xvals, with the column named
Wr.Hnd, which contains only two values of interest-the handspans of 14.5
and 24 cm.
R> xvals <- data.frame(Wr.Hnd=c(14.5,24))
R> xvals
Wr.Hnd
1 14.5
2 24.0
Now, when predict is called, the first argument must be the fitted model
object of interest, survfit for this example. Next, in the argument newdata,
you pass the specially constructed data frame containing the specified predictor values. To the interval argument you must specify "confidence" as
a character string value. The confidence level, here set for 95 percent, is
passed (on the scale of a probability) to level.
R> mypred.ci <- predict(survfit,newdata=xvals,interval="confidence",level=0.95)
R> mypred.ci
fit lwr upr
1 159.1446 156.4956 161.7936
2 188.7524 185.5726 191.9323
This call will return a matrix with three columns, whose number (and
order) of rows correspond to the predictor values you supplied in the newdata
data frame. The first column, with a heading of fit, is the point estimate on
the regression line; you can see that these numbers match the values you
worked out earlier. The other columns provide the lower and upper CI limits as the lwr and upr columns, respectively. In this case, you'd interpret this
as 95 percent confidence that the mean height of a student with a handspan
of 14.5 cm lies somewhere between 156.5 cm and 161.8 cm and lies between
185.6 cm and 191.9 cm for a handspan of 24 cm (when rounded to 1 d.p.).
Remember, these CIs, calculated as per (20.5) through predict, are for the
mean response value.
Prediction Intervals for Individual Observations
The predict function will also provide your prediction intervals. To find the
prediction interval for possible individual observations with a certain probability, you simply need to change the interval argument to "prediction".
R> mypred.pi <- predict(survfit,newdata=xvals,interval="prediction",level=0.95)
R> mypred.pi
Simple Linear Regression 463
fit lwr upr
1 159.1446 143.3286 174.9605
2 188.7524 172.8390 204.6659
Notice that the fitted values remain the same, as Equations (20.5) and
(20.6) indicate. The widths of the PIs, however, are significantly larger than
those of the corresponding CIs-this is because raw observations themselves,
at a specific x value, will naturally be more variable than their mean.
Interpretation changes accordingly. The intervals describe where raw
student heights are predicted to lie "95 percent of the time." For a handspan
of 14.5 cm, the model predicts individual observations to lie somewhere
between 143.3 cm and 175.0 cm with a probability of 0.95; for a handspan
of 24 cm, the same PI is estimated at 172.8 cm and 204.7 cm (when rounded
to 1 d.p.).
20.4.3 Plotting Intervals
Both CIs and PIs are well suited to visualization for simple linear regression
models. With the following code, you can start off Figure 20-3 by plotting the
data and estimated regression line just as for Figure 20-2, but this time using
xlim and ylim in plot to widen the x- and y-limits a little in order to accommodate the full length and breadth of the CI and PI.
R> plot(survey$Height~survey$Wr.Hnd,xlim=c(13,24),ylim=c(140,205),
xlab="Writing handspan (cm)",ylab="Height (cm)")
R> abline(survfit,lwd=2)
To this you add the locations of the fitted values for x = 14.5 and x = 24,
as well as two sets of vertical lines showing the CIs and PIs.
R> points(xvals[,1],mypred.ci[,1],pch=8)
R> segments(x0=c(14.5,24),y0=c(mypred.pi[1,2],mypred.pi[2,2]),
x1=c(14.5,24),y1=c(mypred.pi[1,3],mypred.pi[2,3]),col="gray",lwd=3)
R> segments(x0=c(14.5,24),y0=c(mypred.ci[1,2],mypred.ci[2,2]),
x1=c(14.5,24),y1=c(mypred.ci[1,3],mypred.ci[2,3]),lwd=2)
The call to points marks the fitted values for these two particular values
of x. The first call to segments lays down the PIs as thickened vertical gray
lines, and the second lays down the CIs as the shorter vertical black lines.
The coordinates for these plotted line segments are taken directly from the
mypred.pi and mypred.ci objects, respectively.
You can also produce "bands" around the fitted regression line that
mark one or both of these intervals over all values of the predictor. From
a programming standpoint, this isn't technically possible for a continuous
variable, but you can achieve it practically by defining a fine sequence of values along the x-axis (using seq with a high length value) and evaluating the
CI and PI at every point in this fine sequence. Then you just join resulting
points as lines when plotting.
464 Chapter 20
Figure 20-3: The student height regression example, with a fitted
regression line and point estimates at x = 14.5 and x = 24 and
with corresponding 95 percent CIs (black vertical lines) and PIs
(gray vertical lines). The dashed black and dashed gray lines
provide 95 percent confidence and prediction bands for the
response variable over the visible range of x values.
In R, this requires you to rerun the predict command as follows:
R> xseq <- data.frame(Wr.Hnd=seq(12,25,length=100))
R> ci.band <- predict(survfit,newdata=xseq,interval="confidence",level=0.95)
R> pi.band <- predict(survfit,newdata=xseq,interval="prediction",level=0.95)
The first line in this code creates the fine sequence of predictor values
and stores it in the format required by the newdata argument. The y-axis
coordinates for CI and PI bands are stored as the second and third columns
of the matrix objects ci.band and pi.band. Finally, lines is used to add each of
the four dashed lines corresponding to the upper and lower limits of the two
intervals, and a legend adds a final touch.
R> lines(xseq[,1],ci.band[,2],lty=2)
R> lines(xseq[,1],ci.band[,3],lty=2)
R> lines(xseq[,1],pi.band[,2],lty=2,col="gray")
R> lines(xseq[,1],pi.band[,3],lty=2,col="gray")
R> legend("topleft",legend=c("Fit","95% CI","95% PI"),lty=c(1,2,2),
col=c("black","black","gray"),lwd=c(2,1,1))
Note that the black dashed CI bands meet the vertical black lines and
the gray dashed PI bands meet the vertical gray lines for the two individual x
values from earlier, just as you'd expect.
Simple Linear Regression 465
Figure 20-3 shows the end result of all these additions to the plot. The
"bowing inwards" curvature of the intervals is characteristic of this kind of
plot and is especially visible in the CI. This curve occurs because there is
naturally less variation if you're predicting where there are more data. For
more information on predict for linear model objects, take a look at the
?predict.lm help file.
20.4.4 Interpolation vs. Extrapolation
Before finishing this introduction to prediction, it's important to clarify the definitions of two key terms: interpolation and extrapolation. These
terms describe the nature of a given prediction. A prediction is referred
to as interpolation if the x value you specify falls within the range of your
observed data; extrapolation is when the x value of interest lies outside this
range. From the point-predictions you just made, you can see that the location x = 14.5 is an example of interpolation, and x = 24 is an example of
extrapolation.
In general, interpolation is preferable to extrapolation-it makes more
sense to use a fitted model for prediction in the vicinity of data that have
already been observed. Extrapolation that isn't too far out of that vicinity
may still be considered reliable, though. The extrapolation for the student
height example at x = 24 is a case in point. This is outside the range of the
observed data, but not by much in terms of scale, and the estimated intervals for the expected value of y^ = 188.75 cm appear, at least visually, not
unreasonable given the distribution of the other observations. In contrast,
it would make less sense to use the fitted model to predict student height at
a handspan of, say, 50 cm:
R> predict(survfit,newdata=data.frame(Wr.Hnd=50),interval="confidence",
level=0.95)
fit lwr upr
1 269.7845 251.9583 287.6106
Such an extreme extrapolation suggests that the mean height of an individual with a handspan of 50 cm is almost 270 cm, both being fairly unrealistic measurements. The same is true in the other direction; the intercept ??^
0
doesn't have a particularly useful practical interpretation, indicating that the
mean height of a student with a handspan of 0 cm is around 114 cm.
The main message here is to use common sense when making any prediction from a linear model fit. In terms of the reliability of the results, predictions made at values within an appropriate proximity of the observed data
are preferable.
466 Chapter 20
Exercise 20.1
Continue to use the survey data frame from the package MASS for the
next few exercises.
a. Using your fitted model of student height on writing handspan,
survfit, provide point estimates and 99 percent confidence
intervals for the mean student height for handspans of 12, 15.2,
17, and 19.9 cm.
b. In Section 20.1, you defined the object incomplete.obs, a numeric
vector that provides the records of survey that were automatically removed from consideration when estimating the model
parameters. Now, use the incomplete.obs vector along with survey
and Equation (20.3) to calculate ??^
0 and ??^
1 in R. (Remember
the functions mean, sd, and cor. Ensure your answers match the
output from survfit.)
c. The survey data frame has a number of other variables present
aside from Height and Wr.Hnd. For this exercise, the end aim is
to fit a simple linear model to predict the mean student height,
but this time from their pulse rate, given in Pulse (continue to
assume the conditions listed in Section 20.2 are satisfied).
i. Fit the regression model and produce a scatterplot with the
fitted line superimposed upon the data. Make sure you can
write down the fitted model equation and keep the plot
open.
ii. Identify and interpret the point estimate of the slope, as well
as the outcome of the test associated with the hypotheses
H0 : ??1 = 0; HA : ??1 , 0. Also find a 90 percent CI for the
slope parameter.
iii. Using your model, add lines for 90 percent confidence and
prediction interval bands on the plot from (i) and add a
legend to differentiate between the lines.
iv. Create an incomplete.obs vector for the current "height on
pulse" data. Use that vector to calculate the sample mean of
the height observations that were used for the model fitted
in (i). Then add a perfectly horizontal line to the plot at this
mean (use color or line type options to avoid confusion with
the other lines present). What do you notice? Does the plot
support your conclusions from (ii)?
Simple Linear Regression 467
Next, examine the help file for the mtcars data set, which you first saw
in Exercise 13.4 on page 287. For this exercise, the goal is to model
fuel efficiency, measured in miles per gallon (MPG), in terms of the
overall weight of the vehicle (in thousands of pounds).
d. Plot the data-mpg on the y-axis and wt on the x-axis.
e. Fit the simple linear regression model. Add the fitted line to the
plot from (d).
f. Write down the regression equation and interpret the point
estimate of the slope. Is the effect of wt on mean mpg estimated
to be statistically significant?
g. Produce a point estimate and associated 95 percent PI for a
car that weighs 6,000 lbs. Do you trust the model to predict
observations accurately for this value of the explanatory variable?
Why or why not?
20.5 Understanding Categorical Predictors
So far, you've looked at simple linear regression models that rely on continuous explanatory variables, but it's also possible to use a discrete or categorical explanatory variable, made up of k distinct groups or levels, to model
the mean response. You must be able to make the same assumptions noted
in Section 20.2: that observations are all independent of one another and
residuals are normally distributed with an equal variance. To begin with,
you'll look at the simplest case in which k = 2 (a binary-valued predictor),
which forms the basis of the slightly more complicated situation in which
the categorical predictor has more than two levels (a multilevel predictor: k > 2).
20.5.1 Binary Variables: k = 2
Turn your attention back to Equation (20.1), where the regression model is
specified as Y |X = ??0 + ??1X + o for a response variable Y and predictor X,
and o ??? N(0,??2
). Now, suppose your predictor variable is categorical, with
only two possible levels (binary; k = 2) and observations coded either 0 or 1.
For this case, (20.1) still holds, but the interpretation of the model parameters, ??0 and ??1, isn't really one of an "intercept" and a "slope" anymore.
Instead, it's better to think of them as being something like two intercepts,
where ??0 provides the baseline or reference value of the response when X = 0
and ??1 represents the additive effect on the mean response if X = 1. In other
words, if X = 0, then Y = ??0 + o; if X = 1, then Y = ??0 + ??1 + o. As usual,
estimation is in terms of finding the mean response y^ ??? E[Y |X = x] as per
Equation (20.2), so the equation becomes y^ = ??^
0 + ??^
1 x.
468 Chapter 20
Go back to the survey data frame and note that you have a Sex variable,
where the students recorded their gender. Look at the documentation on
the help page ?survey or enter something like this:
R> class(survey$Sex)
[1] "factor"
R> table(survey$Sex)
Female Male
118 118
You'll see that the sex data column is a factor vector with two levels,
Female and Male, and that there happens to be an equal number of the two
(one of the 237 records has a missing value for this variable).
You're going to determine whether there is statistical evidence that the
height of a student is affected by sex. This means that you're again interested in modeling height as the response variable, but this time, it's with the
categorical sex variable as the predictor.
To visualize the data, if you make a call to plot as follows, you'll get a
pair of boxplots.
R> plot(survey$Height~survey$Sex)
This is because the response variable specified to the left of the ~
is numeric and the explanatory variable to the right is a factor, and the
default behavior of R in that situation is to produce side-by-side boxplots.
To further emphasize the categorical nature of the explanatory variable, you can superimpose the raw height and sex observations on top of
the boxplots. To do this, just convert the factor vector to numeric with a call
to as.numeric; this can be done directly in a call to points.
R> points(survey$Height~as.numeric(survey$Sex),cex=0.5)
Remember that boxplots mark off the median as the central bold line
but that least-squares linear regressions are defined by the mean outcome,
so it's useful to also display the mean heights according to sex.
R> means.sex <- tapply(survey$Height,INDEX=survey$Sex,FUN=mean,na.rm=TRUE)
R> means.sex
Female Male
165.6867 178.8260
R> points(1:2,means.sex,pch=4,cex=3)
You were introduced to tapply in Section 10.2.3; in this call, the argument na.rm=TRUE is matched to the ellipsis in the definition of tapply and is
passed to mean (you need it to ensure the missing values present in the data
Simple Linear Regression 469
do not end up producing NAs as the results). A further call to points adds
those coordinates (as × symbols) to the image; Figure 20-4 gives the final
result.
Figure 20-4: Boxplots of the student heights split by sex,
with the raw observations and sample means (small ??? and
large × symbols, respectively) superimposed
The plot indicates, overall, that males tend to be taller than females-
but is there statistical evidence of a difference to back this up?
Linear Regression Model of Binary Variables
To answer this with a simple linear regression model, you can use lm to produce least-squares estimates just like with every other model you've fitted
so far.
R> survfit2 <- lm(Height~Sex,data=survey)
R> summary(survfit2)
Call:
lm(formula = Height ~ Sex, data = survey)
Residuals:
Min 1Q Median 3Q Max
-23.886 -5.667 1.174 4.358 21.174
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 165.687 0.730 226.98 <2e-16 ***
SexMale 13.139 1.022 12.85 <2e-16 ***
---
470 Chapter 20
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 7.372 on 206 degrees of freedom
(29 observations deleted due to missingness)
Multiple R-squared: 0.4449, Adjusted R-squared: 0.4422
F-statistic: 165.1 on 1 and 206 DF, p-value: < 2.2e-16
However, because the predictor is a factor vector instead of a numeric
vector, the reporting of the coefficients is slightly different. The estimate of
??0 is again reported as (Intercept); this is the estimate of the mean height
if a student is female. The estimate of ??1 is reported as SexMale. The corresponding regression coefficient of 13.139 is the estimated difference that
is imparted upon the mean height of a student if male. If you look at the
corresponding regression equation
y^ = ??^
0 + ??^
1 x = 165.687 + 13.139x (20.7)
you can see that the model has been fitted assuming the variable x is defined
as "the individual is male"-0 for no/false, 1 for yes/true. In other words,
the level of "female" for the sex variable is assumed as a reference, and it
is the effect of "being male" on mean height that is explicitly estimated. The
hypothesis test for ??0 and ??1 is performed with the same hypotheses defined
in Section 20.3.2:
H0 : ??j = 0
HA : ??j , 0
Again, it's the test for ??1 that's generally of the most interest since it's
this value that tells you whether there is statistical evidence that the mean
response variable is affected by the explanatory variable, that is, if ??1 is significantly different from zero.
Predictions from a Binary Categorical Variable
Because there are only two possible values for x, prediction is straightforward here. When you evaluate the equation, the only decision that needs
to be made is whether ??^
1 needs to be used (in other words, if an individual is male) or not (if an individual is female). For example, you can enter
the following code to create a factor of five extra observations with the same
level names as the original data and store the new data in extra.obs:
R> extra.obs <- factor(c("Female","Male","Male","Male","Female"))
R> extra.obs
[1] Female Male Male Male Female
Levels: Female Male
Then, use predict in the now-familiar fashion to find the mean heights
at those extra values of the predictor. (Remember that when you pass in new
data to predict using the newdata argument, the predictors must be in the
same form as the data that were used to fit the model in the first place.)
Simple Linear Regression 471
R> predict(survfit2,newdata=data.frame(Sex=extra.obs),interval="confidence",
level=0.9)
fit lwr upr
1 165.6867 164.4806 166.8928
2 178.8260 177.6429 180.0092
3 178.8260 177.6429 180.0092
4 178.8260 177.6429 180.0092
5 165.6867 164.4806 166.8928
You can see from the output that the predictions are different only
between the two sets of values-the point estimates of the two instances of
Female are identical, simply ??^
0 with 90 percent CIs. The point estimates and
CIs for the instances of Male are also all the same as each other, based on a
point estimate of ??^
0 + ??^
1.
On its own, admittedly, this example isn't too exciting. However, it's
critical to understand how R presents regression results when using categorical predictors, especially when considering multiple regression in
Chapter 21.
20.5.2 Multilevel Variables: k > 2
It's common to work with data where the categorical predictor variables
have more than two levels so that (k > 2). These can also be referred to as
multilevel categorical variables. To deal with this more complicated situation
while retaining interpretability of your parameters, you must first dummy
code your predictor into k ??? 1 binary variables.
Dummy Coding Multilevel Variables
To see how this is done, assume that you want to find the value of response
variable Y when given the value of a categorical variable X, where X has
k > 2 levels (also assume the conditions for validity of the linear regression
model-Section 20.2-are satisfied).
In regression modeling, dummy coding is the procedure used to create
several binary variables from a categorical variable like X. Instead of the
single categorical variable with possible realizations
X = 1,2,3,. . . , k
you recode it into several yes/no variables-one for each level-with
possible realizations:
X(1) = 0,1; X(2) = 0,1; X(3) = 0,1; . . . ; X(k) = 0,1
As you can see, X(i) represents a binary variable for the ith level of the
original X. For example, if an individual has X = 2 in the original categorical variable, then X(2) = 1 (yes) and all of the others (X(1)
, X(3)
, . . . , X(k)) will
be zero (no).
472 Chapter 20
Suppose X is a variable that can take any one of the k = 4 values 1,
2, 3, or 4, and you've made six observations of this variable: 1, 2, 2, 1, 4, 3.
Table 20-1 shows these observations and their dummy-coded equivalents
X(1)
, X(2)
, X(3)
, and X(4)
.
Table 20-1: Illustrative Example of
Dummy Coding for Six Observations of a
Categorical Variable with k = 4 Groups
X X(1) X(2) X(3) X(4)
1 1 0 0 0
2 0 1 0 0
2 0 1 0 0
1 1 0 0 0
4 0 0 0 1
3 0 0 1 0
In fitting the subsequent model, you usually only use k ??? 1 of the dummy
binary variables-one of the variables acts as a reference or baseline level, and
it's incorporated into the overall intercept of the model. In practice, you
would end up with an estimated model like this,
y^ = ??^
0 + ??^
1X(2) + ??^
2X(3) + . . . + ??^
k???1X(k) (20.8)
assuming 1 is the reference level. As you can see, in addition to the overall
intercept term ??^
0, you have k ???1 other estimated intercept terms that modify
the baseline coefficient ??^
0, depending on which of the original categories an
observation takes on. For example, in light of the coding imposed in (20.8),
if an observation has X(3) = 1 and all other binary values are therefore zero
(so that observation would've had a value of X = 3 for the original categorical variable), the predicted mean value of the response would be y^ = ??^
0 + ??^
2.
On the other hand, because the reference level is defined as 1, if an observation has values of zero for all the binary variables, it implies the observation
originally had X = 1, and the prediction would be simply y^ = ??^
0.
The reason it's necessary to dummy code for categorical variables of
this nature is that, in general, categories cannot be related to each other in
the same numeric sense as continuous variables. It's often not appropriate,
for example, to think that an observation in category 4 is "twice as much"
as one in category 2, which is what the estimation methods would assume.
Binary presence/absence variables are valid, however, and can be easily
incorporated into the modeling framework. Choosing the reference level is
generally of secondary importance-the specific values of the estimated coefficients will change accordingly, but any overall interpretations you make
based on the fitted model will be the same regardless.
Simple Linear Regression 473
NOTE Implementing this dummy-coding approach is technically a form of multiple regression since you're now including several binary variables in the model. It's important,
however, to be aware of the somewhat artificial nature of dummy coding-you should
still think of the multiple coefficients as representing a single categorical variable since
the binary variables X(1)
, . . . , X(k) are not independent of one another. This is why
I've chosen to define these models in this chapter; multiple regression will be formally
discussed in Chapter 21.
Linear Regression Model of Multilevel Variables
R makes working with categorical predictors in this way quite simple since
it automatically dummy codes for any such explanatory variable when you
call lm. There are two things you should check before fitting your model,
though.
1. The categorical variable of interest should be stored as a (formally
unordered) factor.
2. You should check that you're happy with the category assigned as the
reference level (for interpretative purposes-see Section 20.5.3).
You must also of course be happy with the validity of the familiar
assumptions of normality and independence of o.
To demonstrate all these definitions and ideas, let's return to the student survey data from the MASS package and keep "student height" as the
response variable of interest. Among the data is the variable Smoke. This
variable describes the kind of smoker each student reports themselves as,
defined by frequency and split into four categories: "heavy," "never," "occasional," and "regular."
R> is.factor(survey$Smoke)
[1] TRUE
R> table(survey$Smoke)
Heavy Never Occas Regul
11 189 19 17
R> levels(survey$Smoke)
[1] "Heavy" "Never" "Occas" "Regul"
Here, the result from is.factor(survey$Smoke) indicates that you do
indeed have a factor vector at hand, the call to table yields the number of
students in each of the four categories, and as per Chapter 5, you can explicitly request the levels attribute of any R factor via levels.
Let's ask whether there's statistical evidence to support a difference in
mean student height according to smoking frequency. You can create a set
of boxplots of these data with the following two lines; Figure 20-5 shows the
result.
R> boxplot(Height~Smoke,data=survey)
R> points(1:4,tapply(survey$Height,survey$Smoke,mean,na.rm=TRUE),pch=4)
474 Chapter 20
Figure 20-5: Boxplots of the observed student heights split by
smoking frequency; respective sample means marked with ×
Note from earlier R output that unless explicitly defined at creation, the
levels of a factor appear in alphabetical order by default-as is the case for
Smoke-and R will automatically set the first one (as shown in the output of a
call to levels) as the reference level when that factor is used as a predictor in
subsequent model fitting. Fitting the linear model in mind using lm, you can
see from a subsequent call to summary that indeed the first level of Smoke, for
"heavy", has been used as the reference:
R> survfit3 <- lm(Height~Smoke,data=survey)
R> summary(survfit3)
Call:
lm(formula = Height ~ Smoke, data = survey)
Residuals:
Min 1Q Median 3Q Max
-25.02 -6.82 -1.64 8.18 28.18
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 173.7720 3.1028 56.005 <2e-16 ***
SmokeNever -1.9520 3.1933 -0.611 0.542
SmokeOccas -0.7433 3.9553 -0.188 0.851
SmokeRegul 3.6451 4.0625 0.897 0.371
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 9.812 on 205 degrees of freedom
(28 observations deleted due to missingness)
Simple Linear Regression 475
Multiple R-squared: 0.02153, Adjusted R-squared: 0.007214
F-statistic: 1.504 on 3 and 205 DF, p-value: 0.2147
As outlined in Equation (20.8), you get estimates of coefficients corresponding to the dummy binary variables for three of the four possible categories in this example-the three nonreference levels. The observation
in the reference category Heavy is represented solely by ??^
0, designated first
as the overall (Intercept), with the other coefficients providing the effects
associated with an observation in one of the other categories.
Predictions from a Multilevel Categorical Variable
You find point estimates through prediction, as usual.
R> one.of.each <- factor(levels(survey$Smoke))
R> one.of.each
[1] Heavy Never Occas Regul
Levels: Heavy Never Occas Regul
R> predict(survfit3,newdata=data.frame(Smoke=one.of.each),
interval="confidence",level=0.95)
fit lwr upr
1 173.7720 167.6545 179.8895
2 171.8200 170.3319 173.3081
3 173.0287 168.1924 177.8651
4 177.4171 172.2469 182.5874
Here, I've created the object one.of.each for illustrative purposes; it
represents one observation in each of the four categories, stored as an
object matching the class (and levels) of the original Smoke data. A student
in the Occas category, for example, is predicted to have a mean height of
173.772 ??? 0.7433 = 173.0287.
The output from the model summary earlier, however, shows that none
of the binary dummy variable coefficients are considered statistically significant from zero (because all the p-values are too large). The results indicate, as you might have suspected, that there's no evidence that smoking
frequency (or more specifically, having a smoking frequency that's different
from the reference level) affects mean student heights based on this sample
of individuals. As is common, the baseline coefficient ??^
0 is highly statistically
significant-but that only suggests that the overall intercept probably isn't
zero. (Because your response variable is a measurement of height and will
clearly not be centered anywhere near 0 cm, that result makes sense.) The
confidence intervals supplied are calculated in the usual t-based fashion.
The small R-Squared value reinforces this conclusion, indicating that
barely any of the variation in the response can be explained by changing
the category of smoking frequency. Furthermore, the overall F -test p-value
is rather large at around 0.215, suggesting an overall nonsignificant effect of
the predictor on the response; you'll look at this in more detail in a moment
in Section 20.5.5 and later on in Section 21.3.5.
476 Chapter 20
As noted earlier, it's important that you interpret these results-indeed
any based on a k-level categorical variable in regression-in a collective fashion. You can claim only that there is no discernible effect of smoking on
height because all the p-values for the binary dummy coefficients are nonsignificant. If one of the levels was in fact highly significant (through a small
p-value), it would imply that the smoking factor as defined here, as a whole,
does have a statistically detectable effect on the response (even if the other
two levels were still associated with very high p-values). This will be discussed
further in several more examples in Chapter 21.
20.5.3 Changing the Reference Level
Sometimes you might decide to change the automatically selected reference
level, compared to which the effects of taking on any of the other levels are
estimated. Changing the baseline will result in the estimation of different
coefficients, meaning that individual p-values are subject to change, but
the overall result (in terms of global significance of the factor) will not be
affected. Because of this, altering the reference level is only done for interpretative purposes-sometimes there's an intuitively natural baseline of the
predictor (for example, "Placebo" versus "Drug A" and "Drug B" as a treatment variable in the analysis of some clinical trial) from which you want to
estimate deviation in the mean response with respect to the other possible
categories.
Redefining the reference level can be achieved quickly using the built-in
relevel function in R. This function allows you to choose which level comes
first in the definition of a given factor vector object and will therefore be
designated as the reference level in subsequent model fitting. In the current
example, let's say you'd rather have the nonsmokers as the reference level.
R> SmokeReordered <- relevel(survey$Smoke,ref="Never")
R> levels(SmokeReordered)
[1] "Never" "Heavy" "Occas" "Regul"
The relevel function has moved the Never category into the first position in the new factor vector. If you go ahead fit the model again using
SmokeReordered instead of the original Smoke column of survey, it'll provide
estimates of coefficients associated with the three different levels of smokers.
It's worth noting the differences in the treatment of unordered versus
ordered factor vectors in regression applications. It might seem sensible
to formally order the smoking variable by, for example, increasing the frequency of smoking when creating a new factor vector. However, when an
ordered factor vector is supplied in a call to lm, R reacts in a different way-it
doesn't perform the relatively simple dummy coding discussed here, where
an effect is associated with each optional level to the baseline (technically
referred to as orthogonal contrasts). Instead, the default behavior is to fit the
model based on something called polynomial contrasts, where the effect of the
ordered categorical variable on the response is defined in a more complicated functional form. That discussion is beyond the scope of this text, but
Simple Linear Regression 477
it suffices to say that this approach can be beneficial when your interest lies
in the specific functional nature of "moving up" through an ordered set of
categories. For more on the technical details, see Kuhn and Johnson (2013).
For all relevant regression examples in this book, we'll work exclusively with
unordered factor vectors.
20.5.4 Treating Categorical Variables as Numeric
The way in which lm decides to define the parameters of the fitted model
depends primarily on the kind of data you pass to the function. As discussed,
lm imposes dummy coding only if the explanatory variable is an unordered
factor vector.
Sometimes the categorical data you want to analyze haven't been stored
as a factor in your data object. If the categorical variable is a character vector, lm will implicitly coerce it into a factor. If, however, the intended categorical variable is numeric, then lm performs linear regression exactly as if it
were a continuous numeric predictor; it estimates a single regression coefficient, which is interpreted as a "per-one-unit-change" in the mean response.
This may seem inappropriate if the original explanatory variable is supposed to be made up of distinct groups. In some settings, however, especially
when the variable can be naturally treated as numeric-discrete, this treatment is not only valid statistically but also helps with interpretation.
Let's take a break from the survey data and go back to the ready-to-use
mtcars data set. Say you're interested in the variables mileage, mpg (continuous), and number of cylinders, cyl (discrete; the data set contains cars with
either 4, 6, or 8 cylinders). Now, it's perfectly sensible to automatically think
of cyl as a categorical variable. Taking mpg to be the response variable, boxplots are well suited to reflect the grouped nature of cyl as a predictor; the
result of the following line is given on the left of Figure 20-6:
R> boxplot(mtcars$mpg~mtcars$cyl,xlab="Cylinders",ylab="MPG")
When fitting the associated regression model, you must be aware of what
you're instructing R to do. Since the cyl column of mtcars is numeric, and
not a factor vector per se, lm will treat it as continuous if you just directly
access the data frame.
R> class(mtcars$cyl)
[1] "numeric"
R> carfit <- lm(mpg~cyl,data=mtcars)
R> summary(carfit)
Call:
lm(formula = mpg ~ cyl, data = mtcars)
Residuals:
Min 1Q Median 3Q Max
-4.9814 -2.1185 0.2217 1.0717 7.5186
478 Chapter 20
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 37.8846 2.0738 18.27 < 2e-16 ***
cyl -2.8758 0.3224 -8.92 6.11e-10 ***
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 3.206 on 30 degrees of freedom
Multiple R-squared: 0.7262, Adjusted R-squared: 0.7171
F-statistic: 79.56 on 1 and 30 DF, p-value: 6.113e-10
Just as in earlier sections, you've received an intercept and a slope estimate; the latter is highly statistically significant, indicating that there is evidence against the true value of the slope being zero. Your fitted regression
line is
y^ = ??^
0 + ??^
1 x = 37.88 ??? 2.88x
where y^ is the average mileage and x is numeric-the number of cylinders. For each single additional cylinder, the model says your mileage will
decrease by 2.88 MPG, on average.
It's important to recognize the fact that you've fitted a continuous line
to what is effectively categorical data. The right panel of Figure 20-6, created
with the following lines, highlights this fact:
R> plot(mtcars$mpg~mtcars$cyl,xlab="Cylinders",ylab="MPG")
R> abline(carfit,lwd=2)
Figure 20-6: Left: Boxplots of mileage split by cylinders for the mtcars data set.
Right: Scatterplot of the same data with fitted regression line (treating cyl as
numeric-continuous) superimposed.
Simple Linear Regression 479
Some researchers fit categorical or discrete predictors as continuous
variables purposefully. First, it allows interpolation; for example, you could
use this model to evaluate the average MPG for a 5-cylinder car. Second, it
means there are fewer parameters that require estimation; in other words,
instead of k ??? 1 intercepts for a categorical variable with k groups, you
need only one parameter for the slope. Finally, it can be a convenient way
to control for so-called nuisance variables; this will become clearer in Chapter 21. On the other hand, it means that you no longer get group-specific
information. It can be misleading to proceed in this way if any differences
in the mean response according to the predictor category of an observation
are not well represented linearly-detection of significant effects can be lost
altogether.
At the very least, it's important to recognize this distinction when fitting
models. If you had only just now recognized that R had fitted the cyl variable as continuous and wanted to actually fit the model with cyl as categorical, you'd have to explicitly convert it into a factor vector beforehand or in
the actual call to lm.
R> carfit <- lm(mpg~factor(cyl),data=mtcars)
R> summary(carfit)
Call:
lm(formula = mpg ~ factor(cyl), data = mtcars)
Residuals:
Min 1Q Median 3Q Max
-5.2636 -1.8357 0.0286 1.3893 7.2364
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 26.6636 0.9718 27.437 < 2e-16 ***
factor(cyl)6 -6.9208 1.5583 -4.441 0.000119 ***
factor(cyl)8 -11.5636 1.2986 -8.905 8.57e-10 ***
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 3.223 on 29 degrees of freedom
Multiple R-squared: 0.7325, Adjusted R-squared: 0.714
F-statistic: 39.7 on 2 and 29 DF, p-value: 4.979e-09
Here, by wrapping cyl in a call to factor when specifying the formula
for lm, you can see you've obtained regression coefficient estimates for the
levels of cyl corresponding to 6- and 8-cylinder cars (with the reference level
automatically set to 4-cylinder cars).
480 Chapter 20
20.5.5 Equivalence with One-Way ANOVA
There's one final observation to make about regression models with a
single nominal categorical predictor. Think about the fact that these models describe a mean response value for the k different groups. Does this
remind you of anything? In this particular setting, you're actually doing the
same thing as in one-way ANOVA (Section 19.1): comparing more than two
means and determining whether there is statistical evidence that at least one
mean is different from the others. You need to be able to make the same key
assumptions of independence and normality for both techniques.
In fact, simple linear regression with a single categorical predictor,
implemented using least-squares estimation, is just another way to perform
one-way ANOVA. Or, perhaps more concisely, ANOVA is a special case of
least-squares regression. The outcome of a one-way ANOVA test is a single
p-value quantifying a level of statistical evidence against the null hypothesis
that states that group means are equal. When you have one categorical predictor in a regression, it's exactly that p-value that's reported at the end of
the summary of an lm object-something I've referred to a couple of times now
as the "overall" or "global" significance test (for example, in Section 20.3.3).
Look back to the final result of that global significance test for the
student height modeled by smoking status example-you had a p-value
of 0.2147. This came from an F test statistic of 1.504 with df1 = 3 and
df2 = 205. Now, suppose you were just handed the data and asked to perform a one-way ANOVA of height on smoking. Using the aov function as
introduced in Section 19.1, you'd call something like this:
R> summary(aov(Height~Smoke,data=survey))
Df Sum Sq Mean Sq F value Pr(>F)
Smoke 3 434 144.78 1.504 0.215
Residuals 205 19736 96.27
28 observations deleted due to missingness
Those same values are returned here; you can also find the square root
of the MSE:
R> sqrt(96.27)
[1] 9.811728
This is in fact the "residual standard error" given in the lm summary.
The two conclusions you'd draw about the impact of smoking status on
height (one for the lm output, the other for the ANOVA test) are of course
also the same.
The global test that lm provides isn't just there for the benefit of confirming ANOVA results. As a generalization of ANOVA, least-squares regression models provide more than just coefficient-specific tests. That global test
Simple Linear Regression 481
is formally referred to as the omnibus F-test, and while it is indeed equivalent
to one-way ANOVA in the "single categorical predictor" setting, it's also a
useful overall, stand-alone test of the statistical contribution of several predictors to the outcome value. You'll explore this further in Section 21.3.5
after you've begun modeling your response variable using multiple explanatory variables.
Exercise 20.2
Continue using the survey data frame from the package MASS for the
next few exercises.
a. The survey data set has a variable named Exer, a factor with k = 3
levels describing the amount of physical exercise time each
student gets: none, some, or frequent. Obtain a count of the
number of students in each category and produce side-by-side
boxplots of student height split by exercise.
b. Assuming independence of the observations and normality as
usual, fit a linear regression model with height as the response
variable and exercise as the explanatory variable (dummy
coding). What's the default reference level of the predictor?
Produce a model summary.
c. Draw a conclusion based on the fitted model from (b)-does it
appear that exercise frequency has any impact on mean height?
What is the nature of the estimated effect?
d. Predict the mean heights of one individual in each of the three
exercise categories, accompanied by 95 percent prediction
intervals.
e. Do you arrive at the same result and interpretation for the
height-by-exercise model if you construct an ANOVA table
using aov?
f. Is there any change to the outcome of (e) if you alter the model
so that the reference level of the exercise variable is "none"?
Would you expect there to be?
Now, turn back to the ready-to-use mtcars data set. One of the variables in this data frame is qsec, described as the time in seconds it
takes to race a quarter mile; another is gear, the number of forward
gears (cars in this data set have either 3, 4, or 5 gears).
g. Using the vectors straight from the data frame, fit a simple linear
regression model with qsec as the response variable and gear as
the explanatory variable and interpret the model summary.
482 Chapter 20
h. Explicitly convert gear to a factor vector and refit the model.
Compare the model summary with that from (g). What do
you find?
i. Explain, with the aid of a relevant plot in the same style as the
right image of Figure 20-6, why you think there is a difference
between the two models (g) and (h).
Important Code in This Chapter
Function/operator Brief description First occurrence
lm Fit linear model Section 20.2.3, p. 455
coef Get estimated coefficients Section 20.2.4, p. 457
summary Summarize linear model Section 20.3.1, p. 458
confint Get CIs for estimated coefficients Section 20.3.2, p. 460
predict Predict from linear model Section 20.4.2, p. 463
relevel Change factor reference level Section 20.5.3, p. 477
Simple Linear Regression 483
21
MULT IPLE L INEAR REGRESS ION
Multiple linear regression is a straightforward generalization of the singlepredictor models discussed in the previous
chapter. It allows you to model your continuous response variable in terms of more than one predictor so you can measure the joint effect of several
explanatory variables on the response variable. In this chapter, you'll see
how to model your response variable in this way, and you'll use R to fit the
model using least-squares. You'll also explore other key statistical aspects of
linear modeling in the R environment, such as transforming variables and
including interactive effects.
Multiple linear regression represents an important part of the practice
of statistics. It lets you control or adjust for multiple sources of influence
on the value of the response, rather than just measuring the effect of one
explanatory variable (in most situations, there is more than one contributor to the outcome measurements). At the heart of this class of methods
is the intention to uncover potentially causal relationships between your
response variable and the (joint) effect of any explanatory variables. In reality, causality itself is extremely difficult to establish, but you can strengthen
any evidence of causality by using a well-designed study supported by sound
data collection and by fitting models that might realistically gauge the relationships present in your data.
21.1 Terminology
Before you look at the theory behind multiple regression models, it's important to have a clear understanding of some terminology associated with
variables.
. A lurking variable influences the response, another predictor, or both,
but goes unmeasured (or is not included) in a predictive model. For
example, say a researcher establishes a link between the volume of trash
thrown out by a household and whether the household owns a trampoline. The potential lurking variable here would be the number of
children in the household-this variable is more likely to be positively
associated with an increase in trash and chances of owning a trampoline. An interpretation that suggests owning a trampoline is a cause of
increased waste would be erroneous.
. The presence of a lurking variable can lead to spurious conclusions
about causal relationships between the response and the other predictors, or it can mask a true cause-and-effect association; this kind of error
is referred to as confounding. To put it another way, you can think of confounding as the entanglement of the effects of one or more predictors
on the response.
. A nuisance or extraneous variable is a predictor of secondary or no interest
that has the potential to confound relationships between other variables
and so affect your estimates of the other regression coefficients. Extraneous variables are included in the modeling as a matter of necessity,
but the specific nature of their influence on the response is not the primary interest of the analysis.
These definitions will become clearer once you begin fitting and interpreting the regression models in Section 21.3. The main message I want to
emphasize here, once more, is that correlation does not imply causation. If
a fitted model finds a statistically significant association between a predictor
(or predictors) and a response, it's important to consider the possibility that
lurking variables are contributing to the results and to attempt to control
any confounding before you draw conclusions. Multiple regression models
allow you to do this.
21.2 Theory
Before you start using R to fit regression models, you'll examine the technical definitions of a linear regression model with multiple predictors.
Here, you'll look at how the models work in a mathematical sense and get
a glimpse of the calculations that happen "behind the scenes" when estimating the model parameters in R.
486 Chapter 21
21.2.1 Extending the Simple Model to a Multiple Model
Rather than having just one predictor, you want to determine the value of
a continuous response variable Y given the values of p > 1 independent
explanatory variables X1, X2, . . ., Xp. The overarching model is defined as
Y = ??0 + ??1X1 + ??2X2 + . . . + ??pXp + o, (21.1)
where ??0,. . . , ??p are the regression coefficients and, as before, you assume
independent, normally distributed residuals o ??? N(0,??) around the mean.
In practice, you have n data records; each record provides values for
each of the predictors Xj
; j = {1, . . ., p}. The model to be fitted is given
in terms of the mean response, conditional upon a particular realization of
the set of explanatory variables
y^ = E[Y |X1 = x1, X2 = x2,. . . , Xp = xp] = ??^
0 + ??^
1 x1 + ??^
2 x2 + . . . + ??^
p xp,
where the ??^
js represent estimates of the regression coefficients.
In simple linear regression, where you have only one predictor variable,
recall that the goal is to find the "line of best fit." The idea of least-squares
estimation for linear models with multiple independent predictors follows
much the same motivation. Now, however, in an abstract sense you can think
of the relationship between response and predictors as a multidimensional
plane or surface. You want to find the surface that best fits your multivariate
data in terms of minimizing the overall squared distance between itself and
the raw response data.
More formally, for your n data records, the ??^
js are found as the values
that minimize the sum
Xn
i=1

yi ??? ( ??0 + ??1 x1,i + ??2 x2,i + . . . + ??^
p xp,i)
	2
, (21.2)
where x j,i
is the observed value of individual i for explanatory variable Xj
and yi
is their response value.
21.2.2 Estimating in Matrix Form
The computations involved in minimizing this squared distance (21.2) are
made much easier by a matrix representation of the data. When dealing with n
multivariate observations, you can write Equation (21.1) as follows,
Y = X · ?? + o,
where Y and o denote n × 1 column matrices such that
Y =
???
???
???
???
???
???
???
???
???
???
???
y1
y2
.
.
.
yn
???
???
???
???
???
???
???
???
???
???
???
and o =
???
???
???
???
???
???
???
???
???
???
???
o1
o2
.
.
.
o n
???
???
???
???
???
???
???
???
???
???
???
.
Multiple Linear Regression 487
Here, yi and oi refer to the response observation and random error
term for the ith individual. The quantity ?? is a (p + 1) × 1 column matrix
of the regression coefficients, and then the observed predictor data for all
individuals and explanatory variables are stored in an n × (p + 1) matrix X,
called the design matrix:
?? =
???
???
???
???
???
???
???
???
???
???
???
??0
??1
.
.
.
??p
???
???
???
???
???
???
???
???
???
???
???
and X =
???
???
???
???
???
???
???
???
???
???
???
1 x1,1 . . . xp,1
1 x1,2 . . . xp,2
.
.
.
.
.
.
.
.
.
.
.
.
1 x1,n . . . xp,n
???
???
???
???
???
???
???
???
???
???
???
The minimization of (21.2) providing the estimated regression coefficient values is then found with the following calculation:
??^ =
???
???
???
???
???
???
???
???
???
???
???
??^
0
??^
1
.
.
.
??^
p
???
???
???
???
???
???
???
???
???
???
???
= (X
???
· X)
???1
· X
???
· Y (21.3)
It's important to note the following:
. The symbol · represents matrix multiplication, the superscript ??? represents the transpose, and ???1
represents the inverse when applied to
matrices (as per Section 3.3).
. Extending the size of ?? and X (note the leading column of 1s in X) to
create structures of size p + 1 (as opposed to just the number of predictors p) allows for the estimation of the overall intercept ??0.
. As well as (21.3), the design matrix plays a crucial role in the estimation
of other quantities, such as the standard errors of the coefficients.
21.2.3 A Basic Example
You can manually estimate the ??j ( j = 0, 1, . . ., p) in R using the functions
covered in Chapter 3: %*% (matrix multiplication), t (matrix transposition),
and solve (matrix inversion). As a quick demonstration, let's say you have
two predictor variables: X1 as continuous and X2 as binary. Your target
regression equation is therefore y^ = ??^
0 + ??^
1 x1 + ??^
2 x2. Suppose you collect
the following data, where the response data, data for X1, and data for X2, for
n = 8 individuals, are given in the columns y, x1, and x2, respectively.
R> demo.data <- data.frame(y=c(1.55,0.42,1.29,0.73,0.76,-1.09,1.41,-0.32),
x1=c(1.13,-0.73,0.12,0.52,-0.54,-1.15,0.20,-1.09),
x2=c(1,0,1,1,0,1,0,1))
R> demo.data
y x1 x2
1 1.55 1.13 1
2 0.42 -0.73 0
3 1.29 0.12 1
488 Chapter 21
4 0.73 0.52 1
5 0.76 -0.54 0
6 -1.09 -1.15 1
7 1.41 0.20 0
8 -0.32 -1.09 1
To get your point estimates in ?? = [??0, ??1, ??2]
??? for the linear model,
you first have to construct X and Y as required by (21.3).
R> Y <- matrix(demo.data$y)
R> Y
[,1]
[1,] 1.55
[2,] 0.42
[3,] 1.29
[4,] 0.73
[5,] 0.76
[6,] -1.09
[7,] 1.41
[8,] -0.32
R> n <- nrow(demo.data)
R> X <- matrix(c(rep(1,n),demo.data$x1,demo.data$x2),nrow=n,ncol=3)
R> X
[,1] [,2] [,3]
[1,] 1 1.13 1
[2,] 1 -0.73 0
[3,] 1 0.12 1
[4,] 1 0.52 1
[5,] 1 -0.54 0
[6,] 1 -1.15 1
[7,] 1 0.20 0
[8,] 1 -1.09 1
Now all you have to do is execute the line corresponding to (21.3).
R> BETA.HAT <- solve(t(X)
R> BETA.HAT
[,1]
[1,] 1.2254572
[2,] 1.0153004
[3,] -0.6980189
You've just used least-squares to fit your model based on the observed
data in demo.data, which results in the estimates ??^
0 = 1.225, ??^
1 = 1.015, and
??^
2 = ???0.698.
Multiple Linear Regression 489
21.3 Implementing in R and Interpreting
Ever helpful, R automatically builds the matrices and carries out all the necessary calculations when you instruct it to fit a multiple linear regression
model. As in simple regression models, you use lm and just include any additional predictors when you specify the formula in the first argument. So that
you can focus on the R syntax and on interpretation, I'll focus only on main
effects for the moment, and then you'll explore more complex relationships
later in the chapter.
When it comes to output and interpretation, working with multiple
explanatory variables follows the same rules as you've seen in Chapter 20.
Any numeric-continuous variables (or a categorical variable being treated
as such) have a slope coefficient that provides a "per-unit-change" quantity.
Any k-group categorical variables (factors, formally unordered) are dummy
coded and provide k ??? 1 intercepts.
21.3.1 Additional Predictors
Let's first confirm the manual matrix calculations from a moment ago.
Using the demo.data object, fit the multiple linear model and examine the
coefficients from that object as follows:
R> demo.fit <- lm(y~x1+x2,data=demo.data)
R> coef(demo.fit)
(Intercept) x1 x2
1.2254572 1.0153004 -0.6980189
You'll see that you obtain exactly the point estimates stored earlier in
BETA.HAT.
With the response variable on the left as usual, you specify the multiple
predictors on the right side of the ~ symbol; altogether this represents the
formula argument. To fit a model with several main effects, use + to separate
any variables you want to include. In fact, you've already seen this notation
in Section 19.2.2, when investigating two-way ANOVA.
To study the interpretation of the parameter estimates of a multiple
linear regression model, let's return to the survey data set in the MASS package. In Chapter 20, you explored several simple linear regression models
based on a response variable of student height, as well as stand-alone predictors of handspan (continuous) and sex (categorical, k = 2). You found that
handspan was highly statistically significant, with the estimated coefficient
suggesting an average increase of about 3.12 cm for each 1 cm increase in
handspan. When you looked at the same t-test using sex as the explanatory
variable, the model also suggested evidence against the null hypothesis, with
"being male" adding around 13.14 cm to the mean height when compared
to the mean for females (the category used as the reference level).
490 Chapter 21
What those models can't tell you is the joint effect of sex and handspan
on predicting height. If you include both predictors in a multiple linear
model, you can (to some extent) reduce any confounding that might otherwise occur in the isolated fits of the effect of either single predictor on
height.
R> survmult <- lm(Height~Wr.Hnd+Sex,data=survey)
R> summary(survmult)
Call:
lm(formula = Height ~ Wr.Hnd + Sex, data = survey)
Residuals:
Min 1Q Median 3Q Max
-17.7479 -4.1830 0.7749 4.6665 21.9253
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 137.6870 5.7131 24.100 < 2e-16 ***
Wr.Hnd 1.5944 0.3229 4.937 1.64e-06 ***
SexMale 9.4898 1.2287 7.724 5.00e-13 ***
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 6.987 on 204 degrees of freedom
(30 observations deleted due to missingness)
Multiple R-squared: 0.5062, Adjusted R-squared: 0.5014
F-statistic: 104.6 on 2 and 204 DF, p-value: < 2.2e-16
The coefficient for handspan is now only about 1.59, almost half of its
corresponding value (3.12 cm) in the stand-alone simple linear regression
for height. Despite this, it's still highly statistically significant in the presence
of sex. The coefficient for sex has also reduced in magnitude when compared with its simple linear model and is also still significant in the presence
of handspan. You'll interpret these new figures in a moment.
As for the rest of the output, the Residual standard error still provides
you with an estimate of the standard error of the random noise term o, and
you're also provided with an R-squared value. When associated with more
than one predictor, the latter is formally referred to as the coefficient of
multiple determination. The calculation of this coefficient, as in the single
predictor setting, comes from the correlations between the variables in the
model. I'll leave the theoretical intricacies to more advanced texts, but it's
important to note that R-squared still represents the proportion of variability
in the response that's explained by the regression; in this example, it sits at
around 0.51.
Multiple Linear Regression 491
You can continue to add explanatory variables in the same way if you
need to do so. In Section 20.5.2, you examined smoking frequency as a
stand-alone categorical predictor for height and found that this explanatory
variable provided no statistical evidence of an impact on the mean response.
But could the smoking variable contribute in a statistically significant way if
you control for handspan and sex?
R> survmult2 <- lm(Height~Wr.Hnd+Sex+Smoke,data=survey)
R> summary(survmult2)
Call:
lm(formula = Height ~ Wr.Hnd + Sex + Smoke, data = survey)
Residuals:
Min 1Q Median 3Q Max
-17.4869 -4.7617 0.7604 4.3691 22.1237
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 137.4056 6.5444 20.996 < 2e-16 ***
Wr.Hnd 1.6042 0.3301 4.860 2.36e-06 ***
SexMale 9.3979 1.2452 7.547 1.51e-12 ***
SmokeNever -0.0442 2.3135 -0.019 0.985
SmokeOccas 1.5267 2.8694 0.532 0.595
SmokeRegul 0.9211 2.9290 0.314 0.753
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 7.023 on 201 degrees of freedom
(30 observations deleted due to missingness)
Multiple R-squared: 0.5085, Adjusted R-squared: 0.4962
F-statistic: 41.59 on 5 and 201 DF, p-value: < 2.2e-16
Since it's a categorical variable with k > 2 levels, Smoke is dummy coded
(with heavy smokers as the default reference level), giving you three extra
intercepts for the three nonreference levels of the variable; the fourth is
incorporated into the overall intercept.
In the summary of the latest fit, you can see that while handspan and sex
continue to yield very small p-values, smoking frequency suggests no such
evidence against the hypotheses of zero coefficients. The smoking variable
has had little effect on the values of the other coefficients compared with the
previous model in survmult, and the R-squared coefficient of multiple determination has barely increased.
One question you might now ask is, if smoking frequency doesn't benefit your ability to predict mean height in any substantial way, should you
remove that variable from the model altogether? This is the primary goal of
model selection: to find the "best" model for predicting the outcome, without
492 Chapter 21
fitting one that is unnecessarily complex (by including more explanatory
variables than is required). You'll look at some common ways researchers
attempt to achieve this in Section 22.2.
21.3.2 Interpreting Marginal Effects
In multiple regression, the estimation of each predictor takes into account
the effect of all other predictors present in the model. A coefficient for a
specific predictor Z should therefore be interpreted as the change in the
mean response for a one-unit increase in Z, while holding all other predictors constant.
As you've determined that smoking frequency still appears to have no
discernible impact on mean height when taking sex and handspan into consideration, return your focus to survmult, the model that includes only the
explanatory variables of sex and handspan. Note the following:
. For students of the same sex (that is, focusing on either just males
or just females), a 1 cm increase in handspan leads to an estimated
increase of 1.5944 cm in mean height.
. For students of similar handspan, males on average will be 9.4898 cm
taller than females.
. The difference in the values of the two estimated predictor coefficients
when compared with their respective simple linear model fits, plus the
fact that both continue to indicate evidence against the null hypothesis of "being zero" in the multivariate fit, suggests that confounding (in
terms of the effect of both handspan and sex on the response variable of
height) is present in the single-predictor models.
The final point highlights the general usefulness of multiple regression.
It shows that, in this example, if you use only single predictor models, the
determination of the "true" impact that each explanatory variable has in
predicting the mean response is misleading since some of the change in
height is determined by sex, but some is also attributed to handspan. It's
worth noting that the coefficient of determination (refer to Section 20.3.3)
for the survmult model is noticeably higher than the same quantity in either
of the single-variate models, so you're actually accounting for more of the
variation in the response by using multiple regression.
The fitted model itself can be thought of as
"Mean height" = 137.687 + 1.594 × "handspan" + 9.49 × "sex" (21.4)
where "handspan" is the writing handspan supplied in centimeters and "sex"
is supplied as either 1 (if male) or 0 (if female).
NOTE The baseline (overall) intercept of around 137.687 cm represents the mean height of a
female with a handspan of 0 cm-again, this is clearly not directly interpretable in the
context of the application. For this kind of situation, some researchers center the offending continuous predictor (or predictors) on zero by subtracting the sample mean of all
the observations on that predictor from each observation prior to fitting the model. The
Multiple Linear Regression 493
centered predictor data are then used in place of the original (untranslated) data. The
resulting fitted model allows you to use the mean value of the untranslated predictor
(in this case handspan) rather than a zero value in order to directly interpret the intercept estimate ??^
0.
21.3.3 Visualizing the Multiple Linear Model
As shown here, "being male" simply changes the overall intercept by around
9.49 cm:
R> survcoefs <- coef(survmult)
R> survcoefs
(Intercept) Wr.Hnd SexMale
137.686951 1.594446 9.489814
R> as.numeric(survcoefs[1]+survcoefs[3])
[1] 147.1768
Because of this, you could also write (21.4) as two equations. Here's the
equation for female students:
"Mean height" = 137.687 + 1.594 × "handspan"
Here's the equation for male students:
"Mean height" = (137.687 + 9.4898) + 1.594 × "handspan"
= 147.177 + 1.594 × "handspan"
This is handy because it allows you to visualize the multivariate model in
much the same way as you can the simple linear models. This code produces
Figure 21-1:
R> plot(survey$Height~survey$Wr.Hnd,
col=c("gray","black")[as.numeric(survey$Sex)],
pch=16,xlab="Writing handspan",ylab="Height")
R> abline(a=survcoefs[1],b=survcoefs[2],col="gray",lwd=2)
R> abline(a=survcoefs[1]+survcoefs[3],b=survcoefs[2],col="black",lwd=2)
R> legend("topleft",legend=levels(survey$Sex),col=c("gray","black"),pch=16)
First, a scatterplot of the height and handspan observations, split by sex,
is drawn. Then, abline adds the line corresponding to females and adds a
second one corresponding to males, based on those two equations.
Although this plot might look like two separate simple linear model
fits, one for each level of sex, it's important to recognize that isn't the case.
You're effectively looking at a representation of a multivariate model on a
two-dimensional canvas, where the statistics that determine the fit of the two
visible lines have been estimated "jointly," in other words, when considering
both predictors.
494 Chapter 21
Figure 21-1: Visualizing the observed data and fitted multiple
linear model of student height modeled by handspan and sex
21.3.4 Finding Confidence Intervals
As in Chapter 20, you can easily find confidence intervals for any of the
regression parameters in multiple regression models with confint. Using
survmult2, the object of the fitted model for student height including the
smoking frequency predictor, the output of a call to confint looks like this:
R> confint(survmult2)
2.5 % 97.5 %
(Intercept) 124.5010442 150.310074
Wr.Hnd 0.9534078 2.255053
SexMale 6.9426040 11.853129
SmokeNever -4.6061148 4.517705
SmokeOccas -4.1312384 7.184710
SmokeRegul -4.8543683 6.696525
Note that the Wr.Hnd and SexMale variables were shown to be statistically
significant at the 5 percent level in the earlier model summary and that their
95 percent confidence levels do not include the null value of zero. On the
other hand, all the coefficients for the dummy variables associated with the
smoking frequency predictor are all nonsignificant, and their confidence
intervals clearly include zero. This reflects the fact that the smoking variable
isn't, as a whole, considered statistically significant in this particular model.
Multiple Linear Regression 495
21.3.5 Omnibus F-Test
First encountered in Section 20.5.2 in the context of multilevel predictors,
you can think of the omnibus F-test more generally for multiple regression
models as a test with the following hypotheses:
H0 : ??1 = ??2 = . . . = ??p = 0
HA : At least one of the ??j , 0 (for j = 1, . . ., p) (21.5)
The test is effectively comparing the amount of error attributed to the
"null" model (in other words, one with an intercept only) with the amount
of error attributed to the predictors when all the predictors are present. In
other words, the more the predictors are able to model the response, the
more error they explain, giving you a more extreme F statistic and therefore
a smaller p-value. The single result makes the test especially useful when you
have many explanatory variables. The test works the same regardless of the
mix of predictors you have in a given model: one or more might be continuous, discrete, binary, and/or categorical with k > 2 levels. When multiple
regression models are fitted, the amount of output alone can take time to
digest and interpret, and care must be taken to avoid Type I errors (incorrect rejection of a true null hypothesis-refer to Section 18.5).
The F-test helps boil all that down, allowing you to conclude either of
the following:
1. Evidence against H0 if the associated p-value is smaller than your chosen
significance level ??, which suggests that your regression-your combination of the explanatory variables-does a significantly better job of
predicting the response than if you removed all those predictors.
2. No evidence against H0 if the associated p-value is larger than ??, which
suggests that using the predictors has no tangible benefit over having an
intercept alone.
The downside is that the test doesn't tell you which of the predictors (or
which subset thereof) is having a beneficial impact on the fit of the model,
nor does it tell you anything about their coefficients or respective standard
errors.
You can compute the F-test statistic using the coefficient of determination, R
2
, from the fitted regression model. Let p be the number of regression parameters requiring estimation, excluding the intercept ??0. Then,
F =
R
2
(n ??? p ??? 1)
(1 ??? R2
)p
, (21.6)
where n is the number of observations used in fitting the model (after
records with missing values have been deleted). Then, under H0 in (21.5),
F follows an F distribution (see Section 16.2.5 and also Section 19.1.2) with
df1 = p, df2 = n???p???1 degrees of freedom. The p-value associated with (21.6)
is yielded as the upper-tail area of that F distribution.
As a quick exercise to confirm this, turn your attention back to the fitted
multiple regression model survmult2 in Section 21.3.1, which is the model
496 Chapter 21
for student height by handspan, sex, and smoking status from survey. You
can extract the coefficient of multiple determination from the summary report
(using the technique noted in Section 20.3.4).
R> R2 <- summary(survmult2)$r.squared
R> R2
[1] 0.508469
This matches the multiple R-squared value from Section 21.3.1. Then,
you can get n as the original size of the data set in survey minus any missing
values (reported as 30 in the earlier summary output).
R> n <- nrow(survey)-30
R> n
[1] 207
You get p as the number of estimated regression parameters (minus 1
for the intercept).
R> p <- length(coef(survmult2))-1
R> p
[1] 5
You can then confirm the value of n ??? p ??? 1, which matches the summary
output (201 degrees of freedom):
R> n-p-1
[1] 201
Finally, you find the test statistic F as dictated by (21.6), and you can
use the pf function as follows to obtain the corresponding p-value for
the test:
R> Fstat <- (R2*(n-p-1))/((1-R2)*p)
R> Fstat
[1] 41.58529
R> 1-pf(Fstat,df1=p,df2=n-p-1)
[1] 0
You can see that the omnibus F-test for this example gives a p-value
that's so small, it's effectively zero. These calculations match the relevant
results reported in the output of summary(survmult2) completely.
Looking back at the student height multiple regression fit based on
handspan, sex, and smoking in survmult2 in Section 21.3.1, it's little surprise
that with two of the predictors yielding small p-values, the omnibus F-test
suggests strong evidence against H0 based on (21.5). This highlights the
"umbrella" nature of the omnibus test: although the smoking frequency
variable itself doesn't appear to contribute anything statistically important,
Multiple Linear Regression 497
the F-test for that model still suggests survmult2 should be preferred over a
"no-predictor" model, because both handspan and sex are important.
21.3.6 Predicting from a Multiple Linear Model
Prediction (or forecasting ) for multiple regression follows the same rules as
for simple regression. It's important to remember that point predictions
found for a particular covariate profile-the collection of predictor values for
a given individual-are associated with the mean (or expected value) of the
response; that confidence intervals provide measures for mean responses;
and that prediction intervals provide measures for raw observations. You also
have to consider the issue of interpolation (predictions based on x values
that fall within the range of the originally observed covariate data) versus
extrapolation (prediction from x values that fall outside the range of said
data). Other than that, the R syntax for predict is identical to that used in
Section 20.4.
As an example, using the model fitted on student height as a linear
function of handspan and sex (in survmult), you can estimate the mean
height of a male student with a writing handspan of 16.5 cm, together with
a confidence interval.
R> predict(survmult,newdata=data.frame(Wr.Hnd=16.5,Sex="Male"),
interval="confidence",level=0.95)
fit lwr upr
1 173.4851 170.9419 176.0283
The result indicates that you have an expected value of about 173.48 cm
and that you can be 95 percent confident the true value lies somewhere
between 170.94 and 176.03 (rounded to 2 d.p.). In the same way, the mean
height of a female with a handspan of 13 cm is estimated at 158.42 cm, with
a 99 percent prediction interval of 139.76 to 177.07.
R> predict(survmult,newdata=data.frame(Wr.Hnd=13,Sex="Female"),
interval="prediction",level=0.99)
fit lwr upr
1 158.4147 139.7611 177.0684
There are in fact two female students in the data set with writing
handspans of 13 cm, as you can see in Figure 21-1. Using your knowledge
of subsetting data frames, you can inspect these two records and select the
three variables of interest.
R> survey[survey$Sex=="Female" & survey$Wr.Hnd==13,c("Sex","Wr.Hnd","Height")]
Sex Wr.Hnd Height
45 Female 13 180.34
152 Female 13 165.00
498 Chapter 21
Now, the second female's height falls well inside the prediction interval,
but the first female's height is significantly higher than the upper limit. It's
important to realize that, technically, nothing has gone wrong here in terms
of the model fitting and interpretation-it's still possible that an observation
can fall outside a prediction interval, even a wide 99 percent interval, though
it's perhaps improbable. There could be any number of reasons for this
occurring. First, the model could be inadequate. For example, you might
be excluding important predictors in the fitted model and therefore have
less predictive power. Second, although the prediction is within the range of
the observed data, it has occurred at one extreme end of the range, where
it's less reliable because your data are relatively sparse. Third, the observation itself may be tainted in some way-perhaps the individual recorded
her handspan incorrectly, in which case her invalid observation should be
removed prior to model fitting. It's with this critical eye that a good statistician will appraise data and models; this is a skill that I'll emphasize further as
this chapter unfolds.
Exercise 21.1
In the MASS package, you'll find the data frame cats, which provides
data on sex, body weight (in kilograms), and heart weight (in grams)
for 144 household cats (see Venables and Ripley, 2002, for further
details); you can read the documentation with a call to ?cats. Load
the MASS package with a call to library("MASS"), and access the object
directly by entering cats at the console prompt.
a. Plot heart weight on the vertical axis and body weight on the
horizontal axis, using different colors or point characters to
distinguish between male and female cats. Annotate your plot
with a legend and appropriate axis labels.
b. Fit a least-squares multiple linear regression model using heart
weight as the response variable and the other two variables as
predictors, and view a model summary.
i. Write down the equation for the fitted model and interpret
the estimated regression coefficients for body weight and
sex. Are both statistically significant? What does this say
about the relationship between the response and predictors?
ii. Report and interpret the coefficient of determination and
the outcome of the omnibus F -test.
c. Tilman's cat, Sigma, is a 3.4 kg female. Use your model to estimate her mean heart weight and provide a 95 percent prediction
interval.
Multiple Linear Regression 499
d. Use predict to superimpose continuous lines based on the fitted
linear model on your plot from (a), one for male cats and one
for female. What do you notice? Does this reflect the statistical
significance (or lack thereof) of the parameter estimates?
The boot package (Davison and Hinkley, 1997; Canty and Ripley,
2015) is another library of R code that's included with the standard
installation but isn't automatically loaded. Load boot with a call to
library("boot"). You'll find a data frame called nuclear, which contains
data on the construction of nuclear power plants in the United States
in the late 1960s (Cox and Snell, 1981).
e. Access the documentation by entering ?nuclear at the prompt
and examine the details of the variables. (Note there is a mistake
for date, which provides the date that the construction permits
were issued-it should read "measured in years since January
1 1900 to the nearest month.") Use pairs to produce a quick
scatterplot matrix of the data.
f. One of the original objectives was to predict the cost of further
construction of these power plants. Create a fit and summary of
a linear regression model that aims to model cost by t1 and t2,
two variables that describe different elapsed times associated with
the application for and issue of various permits. Take note of the
estimated regression coefficients and their significance in the
fitted model.
g. Refit the model, but this time also include an effect for the date
the construction permit was issued. Contrast the output for this
new model against the previous one. What do you notice, and
what does this information suggest about the relationships in the
data with respect to these predictors?
h. Fit a third model for power plant cost, using the predictors for
"date of permit issue," "power plant capacity," and the binary
variable describing whether the plant was sited in the northeastern United States. Write down the fitted model equation
and provide 95 percent confidence intervals for each estimated
coefficient.
The following table gives an excerpt of a historical data set compiled between 1961 and 1973. It concerns the annual murder rate in
Detroit, Michigan; the data were originally presented and analyzed
by Fisher (1976) and are reproduced here from Harraway (1995).
In the data set you'll find the number of murders, police officers,
and gun licenses issued per 100,000 population, as well as the overall
unemployment rate as a percentage of the overall population.
500 Chapter 21
Murders Police Unemployment Guns
8.60 260.35 11.0 178.15
8.90 269.80 7.0 156.41
8.52 272.04 5.2 198.02
8.89 272.96 4.3 222.10
13.07 272.51 3.5 301.92
14.57 261.34 3.2 391.22
21.36 268.89 4.1 665.56
28.03 295.99 3.9 1131.21
31.49 319.87 3.6 837.60
37.39 341.43 7.1 794.90
46.26 356.59 8.4 817.74
47.24 376.69 7.7 583.17
52.33 390.19 6.3 709.59
i. Create your own data frame in your R workspace and produce
a scatterplot matrix. Which of the variables appears to be most
strongly related to the murder rate?
j. Fit a multiple linear regression model using the number of
murders as the response and all other variables as predictors.
Write down the model equation and interpret the coefficients. Is
it reasonable to state that all relationships between the response
and the predictors are causal?
k. Identify the amount of variation in the response attributed to
the joint effect of the three explanatory variables. Then refit the
model excluding the predictor associated with the largest (in
other words, "most nonsignificant") p-value. Compare the new
coefficient of determination with that of the previous model. Is
there much difference?
l. Use your model from (k) to predict the mean number of murders per 100,000 residents, with 300 police officers and 500
issued gun licenses. Compare this to the mean response if there
were no gun licenses issued and provide 99 percent confidence
intervals for both predictions.
21.4 Transforming Numeric Variables
Sometimes, the linear function as strictly defined by the standard regression
equation, (21.1), can be inadequate when it comes to capturing relationships between a response and selected covariates. You might, for example,
observe curvature in a scatterplot between two numeric variables to which
a perfectly straight line isn't necessarily best suited. To a certain extent, the
Multiple Linear Regression 501
requirement that your data exhibit such linear behavior in order for a linear
regression model to be appropriate can be relaxed by simply transforming
(typically in a nonlinear fashion) certain variables before any estimation or
model fitting takes place.
Numeric transformation refers to the application of a mathematical function to your numeric observations in order to rescale them. Finding the
square root of a number and converting a temperature from Fahrenheit to
Celsius are both examples of a numeric transformation. In the context of
regression, transformation is generally applied only to continuous variables
and can be done in any number of ways. In this section, you'll limit your
attention to examples using the two most common approaches: polynomial
and logarithmic transformations. However, note that the appropriateness of
the methods used to transform variables, and any modeling benefits that
might occur, can only really be considered on a case-by-case basis.
Transformation in general doesn't represent a universal solution to solving problems of nonlinearity in the trends in your data, but it can at least
improve how faithfully a linear model is able to represent those trends.
21.4.1 Polynomial
Following on from a comment made earlier, let's say you observe a curved
relationship in your data such that a straight line isn't a sensible choice for
modeling it. In an effort to fit your data more closely, a polynomial or power
transformation can be applied to a specific predictor variable in your regression model. This is a straightforward technique that, by allowing polynomial
curvature in the relationships, allows changes in that predictor to influence
the response in more complex ways than otherwise possible. You achieve
this by including additional terms in the model definition that represent
the impact of progressively higher powers of the variable of interest on the
response.
To clarify the concept of polynomial curvature, consider the following
sequence between ???4 and 4, as well as the simple vectors computed from it:
R> x <- seq(-4,4,length=50)
R> y <- x
R> y2 <- x + x^2
R> y3 <- x + x^2 + x^3
Here, you're taking the original value of x and calculating specific functionals of it. The vector y, as a copy of x, is clearly linear (in technical terms,
this is a "polynomial of order 1"). You assign y2 to take on an additionally
squared valued of x, providing quadratic behavior-a polynomial of order 2.
Lastly, the vector y3 represents the results of a cubic function of the values of
x, with the inclusion of x raised to the power of 3-a polynomial of order 3.
502 Chapter 21
The following three lines of code produce, separately, the plots from left to
right in Figure 21-2.
R> plot(x,y,type="l")
R> plot(x,y2,type="l")
R> plot(x,y3,type="l")
Figure 21-2: Illustrating linear (left), quadratic (middle), and cubic functions (right) of x
Perhaps a bit more generally, let's say you have data for a continuous
predictor, X, that you want to use to model your response, Y. Following
estimation in the usual way, linearly, the simple model is y^ = ??^
0 + ??^
1 x; a
quadratic trend in X can be modeled via the multiple regression y^ = ??^
0 +
??^
1 x+ ??^
2 x
2
; a cubic relationship can be captured by y^ = ??^
0+ ??^
1 x+ ??^
2 x
2+ ??^
3 x
3
;
and so on. From the plots in Figure 21-2, a good way to interpret the effects
of including these extra terms is in the complexity of the curves that can be
captured. At order 1, the linear relationship allows no curvature. At order 2,
a quadratic function of any given variable allows one "bend." At order 3, the
model can cope with two bends in the relationship, and this continues if you
keep adding terms corresponding to increasing powers of the covariate. The
regression coefficients associated with these terms (all implied to be 1 in
the code that produced the previous plots) are able to control the specific
appearance (in other words, the strength and direction) of the curvature.
Fitting a Polynomial Transformation
Return your attention to the built-in mtcars data set. Consider the disp variable, which describes engine displacement volume in cubic inches, against
a response variable of miles per gallon. If you examine a plot of the data in
Figure 21-3, you can see that there does appear to be a slight yet noticeable
curve in the relationship between displacement and mileage.
R> plot(mtcars$disp,mtcars$mpg,xlab="Displacement (cu. in.)",ylab="MPG")
Multiple Linear Regression 503
Figure 21-3: Scatterplot of miles per gallon and engine
displacement, for the mtcars data
Is the straight line that a simple linear regression model would provide
really the best way to represent this relationship? To investigate this, start by
fitting that basic linear setup.
R> car.order1 <- lm(mpg~disp,data=mtcars)
R> summary(car.order1)
Call:
lm(formula = mpg ~ disp, data = mtcars)
Residuals:
Min 1Q Median 3Q Max
-4.8922 -2.2022 -0.9631 1.6272 7.2305
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 29.599855 1.229720 24.070 < 2e-16 ***
disp -0.041215 0.004712 -8.747 9.38e-10 ***
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 3.251 on 30 degrees of freedom
Multiple R-squared: 0.7183, Adjusted R-squared: 0.709
F-statistic: 76.51 on 1 and 30 DF, p-value: 9.38e-10
504 Chapter 21
This clearly indicates statistical evidence of a negative linear impact of
displacement on mileage-for each additional cubic inch of displacement,
the mean response decreases by about 0.041 miles per gallon.
Now, try to capture the apparent curve in the data by adding a quadratic term in disp to the model. You can do this in two ways. First, you could
create a new vector in the workspace by simply squaring the mtcars$disp vector and then supplying the result to the formula in lm. Second, you could
specify disp^2 directly as an additive term in the formula. If you do it this
way, it's essential to wrap that particular expression in a call to I as follows:
R> car.order2 <- lm(mpg~disp+I(disp^2),data=mtcars)
R> summary(car.order2)
Call:
lm(formula = mpg ~ disp + I(disp^2), data = mtcars)
Residuals:
Min 1Q Median 3Q Max
-3.9112 -1.5269 -0.3124 1.3489 5.3946
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 3.583e+01 2.209e+00 16.221 4.39e-16 ***
disp -1.053e-01 2.028e-02 -5.192 1.49e-05 ***
I(disp^2) 1.255e-04 3.891e-05 3.226 0.0031 **
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 2.837 on 29 degrees of freedom
Multiple R-squared: 0.7927, Adjusted R-squared: 0.7784
F-statistic: 55.46 on 2 and 29 DF, p-value: 1.229e-10
Use of the I function around a given term in the formula is necessary
when said term requires an arithmetic calculation-in this case, disp^2-
before the model itself is actually fitted.
Turning to the fitted multiple regression model itself, you can see that
the contribution of the squared component is statistically significant-the
output corresponding to I(disp^2) shows a p-value of 0.0031. This implies
that even if a linear trend is taken into account, the model that includes a
quadratic component (which introduces a curve) is a better-fitting model.
This conclusion is supported by a noticeably higher coefficient of determination compared to the first fit (0.7927 against 0.7183). You can see the fit of
this quadratic curve in Figure 21-4 (code for which follows shortly).
Multiple Linear Regression 505
Here you might reasonably wonder whether you can improve the ability of the model to capture the relationship further by adding yet another
higher-order term in the covariate of interest. To that end:
R> car.order3 <- lm(mpg~disp+I(disp^2)+I(disp^3),data=mtcars)
R> summary(car.order3)
Call:
lm(formula = mpg ~ disp + I(disp^2) + I(disp^3), data = mtcars)
Residuals:
Min 1Q Median 3Q Max
-3.0896 -1.5653 -0.3619 1.4368 4.7617
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 5.070e+01 3.809e+00 13.310 1.25e-13 ***
disp -3.372e-01 5.526e-02 -6.102 1.39e-06 ***
I(disp^2) 1.109e-03 2.265e-04 4.897 3.68e-05 ***
I(disp^3) -1.217e-06 2.776e-07 -4.382 0.00015 ***
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 2.224 on 28 degrees of freedom
Multiple R-squared: 0.8771, Adjusted R-squared: 0.8639
F-statistic: 66.58 on 3 and 28 DF, p-value: 7.347e-13
The output shows that a cubic component also offers a statistically significant contribution. However, if you were to continue adding higher-order
terms, you'd find that fitting a polynomial of order 4 to these data isn't able
to improve the fit at all, with several coefficients being rendered nonsignificant (the order 4 fit isn't shown).
So, letting y^ be miles per gallon and x be displacement in cubic inches,
and expanding the e-notation from the previous output, the fitted multiple
regression model is
y^ = 50.7 ??? 0.3372x + 0.0011x
2 ??? 0.000001x
3
,
which is precisely what the order 3 line in the left panel of Figure 21-4
reflects.
Plotting the Polynomial Fit
To address the plot itself, you visualize the data and the first (simple linear)
model in car.order1 in the usual way. To begin Figure 21-4, execute the following code:
R> plot(mtcars$disp,mtcars$mpg,xlab="Displacement (cu. in.)",ylab="MPG")
R> abline(car.order1)
506 Chapter 21
Figure 21-4: Three different models, polynomials of orders 1, 2, and 3, fitted to the
"mileage per displacement" relationship from the mtcars data set. Left: Visible plot
limits constrained to the data. Right: Visible plot limits widened considerably to
illustrate unreliability in extrapolation.
It's a little more difficult to add the line corresponding to either of the
polynomial-termed models since abline is equipped to handle only straightline trends. One way to do this is to make use of predict for each value in
a sequence that represents the desired values of the explanatory variable.
(I favor this approach because it also allows you to simultaneously calculate confidence and prediction bands if you want.) To add the line for the
order 2 model only, first create the required sequence over the observed
range of disp.
R> disp.seq <- seq(min(mtcars$disp)-50,max(mtcars$disp)+50,length=30)
Here, the sequence has been widened a little by minus and plus 50 to
predict a small amount on either side of the scope of the original covariate
data, so the curve meets the edges of the graph. Then you make the prediction itself and superimpose the fitted line.
R> car.order2.pred <- predict(car.order2,newdata=data.frame(disp=disp.seq))
R> lines(disp.seq,car.order2.pred,lty=2)
You use the same technique, followed by the final addition of the legend, for the order 3 polynomial.
R> car.order3.pred <- predict(car.order3,newdata=data.frame(disp=disp.seq))
R> lines(disp.seq,car.order3.pred,lty=3)
R> legend("topright",lty=1:3,
legend=c("order 1 (linear)","order 2 (quadratic)","order 3 (cubic)"))
The result of all this is on the left panel of Figure 21-4. Even though
you've used raw data from only one covariate, disp, the example illustrated
here is considered multiple regression because more than one parameter
Multiple Linear Regression 507
(in addition to the universal intercept ??0) required estimation in the order
2 and 3 models.
The different types of trend lines fitted to the mileage and displacement
data clearly show different interpretations of the relationship. Visually, you
could reasonably argue that the simple linear fit is inadequate at modeling
the relationship between response and predictor, but it's harder to come to
a clear conclusion when choosing between the order 2 and order 3 versions.
The order 2 fit captures the curve that tapers off as disp increases; the order
3 fit additionally allows for a bump (in technical terms a saddle or inflection),
followed by a steeper downward trend in the same domain.
So, which model is "best"? In this case, the statistical significance of the
parameters suggests that the order 3 model should be preferred. Having
said that, there are other things to consider when choosing between different models, which you'll think about more carefully in Section 22.2.
Pitfalls of Polynomials
One particular drawback associated with polynomial terms in linear regression models is the instability of the fitted trend when trying to perform any
kind of extrapolation. The right plot in Figure 21-4 shows the same three
fitted models (MPG by displacement), but this time with a much wider scale
for displacement. As you can see, the validity of these models is questionable. Though the order 2 and 3 models fit MPG acceptably within the range
of the observed data, if you move even slightly outside the maximum threshold of observed displacement values, the predictions of the mean mileage
go wildly off course. The order 2 model in particular becomes completely
nonsensical, suggesting a rapid improvement in MPG once the engine displacement rises over 500 cubic inches. You must keep this natural mathematical behavior of polynomial functions in mind if you're considering
using higher-order terms in your regression models.
To create this plot, the same code that created the left plot can be used;
you simply use xlim to widen the x-axis range and define the disp.seq object
to a correspondingly wider sequence (in this case, I just set xlim=c(10,1000)
with matching from and to limits in the creation of disp.seq).
NOTE Models like this are still referred to as linear regression models, which might seem a
bit confusing since the fitted trends for higher-order polynomials are clearly nonlinear.
This is because linear regression refers to the fact that the function defining the
mean response is linear in terms of the regression parameters ??0, ??1, . . ., ??p. As such,
any transformation applied to individual variables doesn't affect the linearity of the
function with respect to the coefficients themselves.
21.4.2 Logarithmic
In statistical modeling situations where you have positive numeric observations, it's common to perform a log transformation of the data to
dramatically reduce the overall range of the data and bring extreme
508 Chapter 21
observations closer to a measure of centrality. In that sense, transforming to
a logarithmic scale can help reduce the severity of heavily skewed data (see
Section 15.2.4). In the context of regression modeling, log transformations
can be used to capture trends where apparent curves "flatten off," without
the same kind of instability outside the range of the observed data that you
saw with some of the polynomials.
If you need to refresh your memory on logarithms, turn back to Section 2.1.2; it suffices here to note that the logarithm is the power to which
you must raise a base value in order to obtain an x value. For example, in
3
5 = 243, the logarithm is 5 and 3 is the base, expressed as log3
243 = 5.
Because of the ubiquity of the exponential function in common probability
distributions, statisticians almost exclusively work with the natural log (logarithm to the base e). From here, assume all mentions of the log transformation refer to the natural log.
To briefly illustrate the typical behavior of the log transformation, take a
look at Figure 21-5, achieved with the following:
R> plot(1:1000,log(1:1000),type="l",xlab="x",ylab="",ylim=c(-8,8))
R> lines(1:1000,-log(1:1000),lty=2)
R> legend("topleft",legend=c("log(x)","-log(x)"),lty=c(1,2))
This plots the log of the integers 1 to 1000 against the raw values, as
well as plotting the negative log. You can see the way in which the logtransformed values taper off and flatten out as the raw values increase.
Figure 21-5: The log function applied
to integers 1 to 1000
Fitting the Log Transformation
As noted, one use of the log transformation in regression is to allow this
kind of curvature in situations when a perfectly straight line doesn't suit
the observed relationship. For an illustration, return to the mtcars examples
and consider mileage as a function of both horsepower and transmission
Multiple Linear Regression 509
type (variables hp and am, respectively). Create a scatterplot of MPG against
horsepower, with different colors distinguishing between automatic and
manual cars.
R> plot(mtcars$hp,mtcars$mpg,pch=19,col=c("black","gray")[factor(mtcars$am)],
xlab="Horsepower",ylab="MPG")
R> legend("topright",legend=c("auto","man"),col=c("black","gray"),pch=19)
The plotted points shown in Figure 21-6 suggest that curved trends in
horsepower may be more appropriate than straight-line relationships. Note
that you have to explicitly coerce the binary numeric mtcars$am vector to a
factor here in order to use it as a selector for the vector of two colors. You'll
add the lines in after fitting the linear model.
Figure 21-6: Scatterplot of MPG on horsepower, split by
transmission type, with lines corresponding to a multiple linear
regression using a log-scaled effect of horsepower superimposed
Let's do so using the log transformation of horsepower to try to capture
the curved relationship. Since, in this example, you also want to account
for the potential of transmission type to affect the response, this is included
as an additional predictor variable as usual.
R> car.log <- lm(mpg~log(hp)+am,data=mtcars)
R> summary(car.log)
Call:
lm(formula = mpg ~ log(hp) + am, data = mtcars)
510 Chapter 21
Residuals:
Min 1Q Median 3Q Max
-3.9084 -1.7692 -0.1432 1.4032 6.3865
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 63.4842 5.2697 12.047 8.24e-13 ***
log(hp) -9.2383 1.0439 -8.850 9.78e-10 ***
am 4.2025 0.9942 4.227 0.000215 ***
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 2.592 on 29 degrees of freedom
Multiple R-squared: 0.827, Adjusted R-squared: 0.8151
F-statistic: 69.31 on 2 and 29 DF, p-value: 8.949e-12
The output indicates jointly statistically significant effects of both loghorsepower and transmission type on mileage. Keeping transmission constant, the mean MPG drops by around 9.24 for each additional unit of loghorsepower. Having a manual transmission increases the mean MPG by
roughly 4.2 (estimated in this order owing to the coding of am-0 for automatic, 1 for manual; see ?mtcars). The coefficient of determination shows
82.7 percent of the variation in the response is explained by this regression,
suggesting a satisfactory fit.
Plotting the Log Transformation Fit
To visualize the fitted model, you first need to calculate the fitted values for
all desired predictor values. The following code creates a sequence of horsepower values (minus and plus 20 horsepower) and performs the required
prediction for both transmission types.
R> hp.seq <- seq(min(mtcars$hp)-20,max(mtcars$hp)+20,length=30)
R> n <- length(hp.seq)
R> car.log.pred <- predict(car.log,newdata=data.frame(hp=rep(hp.seq,2),
am=rep(c(0,1),each=n)))
In the above code, since you want to plot predictions for both possible
values of am, when using newdata you need to replicate hp.seq twice. Then,
when you provide values for am to newdata, one series of hp.seq is paired with
an appropriately replicated am value of 0, the other with 1. The result of this
is a vector of predictions of length twice that of hp.seq, car.log.pred, with the
first n elements corresponding to automatic cars and the latter n to manuals.
Now you can add these lines to Figure 21-6 with the following:
R> lines(hp.seq,car.log.pred[1:n])
R> lines(hp.seq,car.log.pred[(n+1):(2*n)],col="gray")
Multiple Linear Regression 511
By examining the scatterplot, you can see that the fitted model appears
to do a good job of estimating the joint relationship between horsepower/
transmission and MPG. The statistical significance of transmission type in
this model directly affects the difference between the two added lines. If
am weren't significant, the lines would be closer together; in that case, the
model would be suggesting that one curve would be sufficient to capture
the relationship. As usual, extrapolation too far outside the range of the
observed predictor data isn't a great idea, though it's less unstable for logtransformed trends than for polynomial functions.
21.4.3 Other Transformations
Transformation can involve more than one variable of the data set and isn't
limited to just predictor variables either. In their original investigation into
the mtcars data, Henderson and Velleman (1981) also noted the presence
of the same curved relationships you've uncovered between the response
and variables such as horsepower and displacement. They argued that it's
preferable to use gallons per mile (GPM) instead of MPG as the response
variable to improve linearity. This would involve modeling a transformation
of MPG, namely, that GPM = 1/MPG.
The authors also commented on the limited influence that both horsepower and displacement have on GPM if the weight of the car is included in
a fitted model, because of the relatively high correlations present among
these three predictors (known as multicollinearity). To address this, the
authors created a new predictor variable calculated as horsepower divided
by weight. This measures, in their words, how "overpowered" a car is-and
they proceeded to use that new predictor instead of horsepower or displacement alone. This is just some of the experimentation that took place in the
search for an appropriate way to model these data.
To this end, however you choose to model your own data, the objective of transforming numeric variables should always be to fit a valid model
that represents the data and the relationships more realistically and accurately. When reaching for this goal, there's plenty of freedom in how you
can transform numeric observations in applications of regression methods.
For a further discussion on transformations in linear regression, Chapter 7
of Faraway (2005) provides an informative introduction.
Exercise 21.2
The following table presents data collected in one of Galileo's
famous "ball" experiments, in which he rolled a ball down a ramp
of different heights and measured how far it traveled from the
base of the ramp. For more on this and other interesting examples,
look at "Teaching Statistics with Data of Historic Significance" by
Dickey and Arnold (1995).
512 Chapter 21
Initial height Distance
1000 573
800 534
600 495
450 451
300 395
200 337
100 253
a. Create a data frame in R based on this table and plot the data
points with distance on the y-axis.
b. Galileo believed there was a quadratic relationship between
initial height and the distance traveled.
i. Fit an order 2 polynomial in height, with distance as the
response.
ii. Fit a cubic (order 3) and a quartic (order 4) model for
these data. What do they tell you about the nature of the
relationship?
c. Based on your models from (b), choose the one that you think
best represents the data and plot the fitted line on the raw data.
Add 90 percent confidence bands for mean distance traveled to
the plot.
The contributed R package faraway contains a large number of data
sets that accompany a linear regression textbook by Faraway (2005).
Install the package and then call library("faraway") to load it. One of
the data sets is trees, which provides data on the dimensions of felled
trees of a certain type (see, for example, Atkinson, 1985).
d. Access the data object at the prompt and plot volume against
girth (the latter along the x-axis).
e. Fit two models with Volume as the response: one quadratic model
in Girth and the other based on log transformations of both
Volume and Girth. Write down the model equations for each and
comment on the similarity (or difference) of the fits in terms of
the coefficient of determination and the omnibus F-test.
f. Use predict to add lines to the plot from (d) for each of the two
models from (e). Use different line types; add a corresponding
legend. Also include 95 percent prediction intervals, with line
types matching those of the fitted values (note that for the model
that involves log transformation of the response and the predictor, any returned values from predict will themselves be on the
log scale; you have to back-transform these to the original scale
Multiple Linear Regression 513
using exp before the lines for that model can be superimposed).
Comment on the respective fits and their estimated prediction
intervals.
Lastly, turn your attention back to the mtcars data frame.
g. Fit and summarize a multiple linear regression model to determine mean MPG from horsepower, weight, and displacement.
h. In the spirit of Henderson and Velleman (1981), use I to refit
the model in (g) in terms of GPM = 1/MPG. Which model
explains a greater amount of variation in the response?
21.5 Interactive Terms
So far, you've looked only at the joint main effects of how predictors affect
the outcome variable (and one-to-one transformations thereof). Now you'll
look at interactions between covariates. An interactive effect between predictors is an additional change to the response that occurs at particular combinations of the predictors. In other words, an interactive effect is present if,
for a given covariate profile, the values of the predictors are such that they
produce an effect that augments the stand-alone main effects associated with
those predictors.
21.5.1 Concept and Motivation
Diagrams such as those found in Figure 21-7 are often used to help explain
the concept of interactive effects. These diagrams show your mean response
value, y^, on the vertical axis, as usual, and a predictor value for the variable
x1 on the horizontal axis. They also show a binary categorical variable x2,
which can be either zero or one. These hypothetical variables are labeled as
such in the images.
Figure 21-7: Concept of an interactive effect between two predictors x1 and x2,
on the mean response value y^. Left: Only main effects of x1 and x2 influence y^.
Right: An interaction between x1 and x2 is needed in addition to their main
effects in order to model y^.
514 Chapter 21
The left diagram shows the limit of the models you've considered so far
in this chapter-that both x1 and x2 affect y^ independently of each other.
The right diagram, however, clearly shows that the effect of x1 on y^ changes
completely depending on the value of x2. On the left, only main effects of x1
and x2 are needed to determine y^; on the right, main effects and an interactive effect between x1 and x2 are present.
NOTE When estimating regression models, you always have to accompany interactions with
the main effects of the relevant predictors, for reasons of interpretability. Since interactions are themselves best understood as an augmentation of the main effects, it
makes no sense to remove the latter and leave in the former.
For a good example of an interaction, think about pharmacology. Interactive effects between medicines are relatively common, which is why health
care professionals often ask about other medicines you might be taking.
Consider statins-drugs commonly used to reduce cholesterol. Users of
statins are told to avoid grapefruit juice because it contains natural chemical compounds that inhibit the efficacy of the enzyme responsible for the
correct metabolization of the drug. If an individual is taking statins and not
consuming grapefruit, you would expect a negative relationship between
cholesterol level and statin use (think about "statin use" either as a continuous or as a categorical dosage variable)-as statin use increases or is affirmative, the cholesterol level decreases. On the other hand, for an individual on
statins who is consuming grapefruit, the nature of the relationship between
cholesterol level and statin use could easily be different-weakened negative, neutral, or even positive. If so, since the effect of the statins on cholesterol changes according to the value of another variable-whether or not
grapefruit is consumed-this would be considered an interaction between
those two predictors.
Interactions can occur between categorical variables, numeric variables, or both. It's most common to find two-way interactions-interactions
between exactly two predictors-which is what you'll focus on in Sections 21.5.2 to 21.5.4. Three-way and higher-order interactive effects are
technically possible but less common, partly because they are difficult to
interpret in a real-world context. You'll consider an example of these in
Section 21.5.5.
21.5.2 One Categorical, One Continuous
Generally, a two-way interaction between a categorical and a continuous
predictor should be understood as effecting a change in the slope of the
continuous predictor with respect to the nonreference levels of the categorical predictor. In the presence of a term for the continuous variable, a
categorical variable with k levels will have k ??? 1 main effect terms, so there
will be a further k ??? 1 interactive terms between all the alternative levels of
the categorical variable and the continuous variable.
Multiple Linear Regression 515
The different slopes for x1 by category of x2 for y^ can be seen clearly on
the right of Figure 21-7. In such a situation, in addition to the main effects
for x1 and x2, there would be one interactive term in the fitted model corresponding to x2 = 1. This defines the additive term needed to change the
slope in x1 for x2 = 0 to the new slope in x1 for x2 = 1.
For an example, let's access a new data set. In Exercise 21.2, you looked
at the faraway package (Faraway, 2005) to access the trees data. In this package, you'll also find the diabetes object-a cardiovascular disease data set
detailing characteristics of 403 African Americans (originally investigated
and reported in Schorling et al., 1997; Willems et al., 1997). Install faraway if
you haven't already and load it with library("faraway"). Restrict your attention to the total cholesterol level (chol-continuous), age of the individual (age-continuous), and body frame type (frame-categorical with k = 3
levels: "small" as the reference level, "medium", and "large"). You can see the
data in Figure 21-8, which will be created momentarily.
You'll look at modeling total cholesterol by age and body frame. It
seems logical to expect that cholesterol is related to both age and body
type, so it makes sense to also consider the possibility that the effect of age
on cholesterol is different for individuals of different body frames. To investigate, let's fit the multiple linear regression and include a two-way interaction between the two variables. In the call to lm, you specify the main effects
first, using + as usual, and then specify an interactive effect of two predictors
by using a colon (:) between them.
R> dia.fit <- lm(chol~age+frame+age:frame,data=diabetes)
R> summary(dia.fit)
Call:
lm(formula = chol ~ age + frame + age:frame, data = diabetes)
Residuals:
Min 1Q Median 3Q Max
-131.90 -26.24 -5.33 22.17 226.11
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 155.9636 12.0697 12.922 < 2e-16 ***
age 0.9852 0.2687 3.667 0.00028 ***
framemedium 28.6051 15.5503 1.840 0.06661 .
framelarge 44.9474 18.9842 2.368 0.01840 *
age:framemedium -0.3514 0.3370 -1.043 0.29768
age:framelarge -0.8511 0.3779 -2.252 0.02490 *
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 42.34 on 384 degrees of freedom
(13 observations deleted due to missingness)
516 Chapter 21
Multiple R-squared: 0.07891, Adjusted R-squared: 0.06692
F-statistic: 6.58 on 5 and 384 DF, p-value: 6.849e-06
Inspecting the estimated model parameters in the output, you can see
a main effect coefficient for age, main effect coefficients for the two levels of
frame (that aren't the reference level), and two further terms for the interactive effect of age with those same nonreference levels.
NOTE There's actually a shortcut to doing this in R-the cross-factor notation. The same
model shown previously could have been fitted by using chol~age*frame in lm; the symbol * between two variables in a formula should be interpreted as "include an intercept,
all main effects, and the interaction." I'll use this shortcut from now on.
The output shows the significance of age and some evidence to support
the presence of a main effect of frame. There's also slight indication of significance of the interaction, though it's weak. Assessing significance in this
case, where one predictor is categorical with k > 2 levels, follows the same
rule as noted in the discussion of multilevel variables in Section 20.5.2-if at
least one of the coefficients is significant, the entire effect should be deemed
significant.
The general equation for the fitted model can be written down directly
from the output.
"Mean total cholesterol" = 155.9636 + 0.9852 × "age"
+ 28.6051 × "medium frame"
+ 44.9474 × "large frame"
??? 0.3514 × "age : medium frame"
??? 0.8511 × "age : large frame" (21.7)
I've used a colon (:) to denote the interactive terms to mirror the R
output.
For the reference level of the categorical predictor, body type "small,"
the fitted model can be written down straight from the output.
"Mean total cholesterol" = 155.9636 + 0.9852 × "age"
For a model with the main effects only, changing body type to "medium"
or "large" would affect only the intercept-you know from Section 20.5 that
the relevant effect is simply added to the outcome. The presence of the
interaction, however, means that in addition to the change in the intercept,
the main effect slope of age must now also be changed according to the relevant interactive term. For an individual with a "medium" frame, the model is
"Mean total cholesterol" = 155.9636 + 0.9852 × "age" + 28.6051
??? 0.3514 × "age"
= 184.5687 + (0.9852 ??? 0.3514) × "age"
= 184.5687 + 0.6338 × "age"
Multiple Linear Regression 517
and for an individual with a "large" frame, the model is
"Mean total cholesterol" = 155.9636 + 0.9852 × "age" + 44.9474
??? 0.8511 × "age"
= 200.911 + (0.9852 ??? 0.8511) × "age"
= 200.911 + 0.1341 × "age"
You can easily calculate these in R by accessing the coefficients of the
fitted model object:
R> dia.coef <- coef(dia.fit)
R> dia.coef
(Intercept) age framemedium framelarge
155.9635868 0.9852028 28.6051035 44.9474105
age:framemedium age:framelarge
-0.3513906 -0.8510549
Next, let's sum the relevant components of this vector. Once you have
the sums, you'll be able to plot the fitted model.
R> dia.small <- c(dia.coef[1],dia.coef[2])
R> dia.small
(Intercept) age
155.9635868 0.9852028
R> dia.medium <- c(dia.coef[1]+dia.coef[3],dia.coef[2]+dia.coef[5])
R> dia.medium
(Intercept) age
184.5686904 0.6338122
R> dia.large <- c(dia.coef[1]+dia.coef[4],dia.coef[2]+dia.coef[6])
R> dia.large
(Intercept) age
200.9109973 0.1341479
The three lines are stored as numeric vectors of length 2, with the intercept first and the slope second. This is the form required by the optional
coef argument of abline, which allows you to superimpose these straight lines
on a plot of the raw data. The following code produces Figure 21-8.
R> cols <- c("black","darkgray","lightgray")
R> plot(diabetes$chol~diabetes$age,col=cols[diabetes$frame],
cex=0.5,xlab="age",ylab="cholesterol")
R> abline(coef=dia.small,lwd=2)
R> abline(coef=dia.medium,lwd=2,col="darkgray")
R> abline(coef=dia.large,lwd=2,col="lightgray")
R> legend("topright",legend=c("small frame","medium frame","large frame"),
lty=1,lwd=2,col=cols)
518 Chapter 21
Figure 21-8: Fitted linear model, main effects, and interaction
for mean total cholesterol by age and body frame
If you examine the fitted model in Figure 21-8, it's clear that inclusion
of an interaction between age and body frame has allowed more flexibility
in the way mean total cholesterol relates to the two predictors. The nonparallel nature of the three plotted lines reflects the concept illustrated in
Figure 21-7.
I walked through this to illustrate how the concept works, but in practice
you don't need to go through all of these steps to find the point estimates
(and any associated confidence intervals). You can predict from a fitted linear model with interactions in the same way as for main-effect-only models
through the use of predict.
21.5.3 Two Categorical
You met the concept of interactions between two categorical explanatory
variables in the introduction to two-way ANOVA in Section 19.2. There, you
uncovered evidence of an interactive effect of wool type and tension on the
mean number of warp breaks in lengths of yarn (based on the ready-to-use
warpbreaks data frame). You then visualized the interaction with an interaction plot (Figure 19-2 on page 447), not unlike the diagrams in Figure 21-7.
Let's implement the same model as the last warpbreaks example in Section 19.2.2 in an explicit linear regression format.
R> warp.fit <- lm(breaks~wool*tension,data=warpbreaks)
R> summary(warp.fit)
Call:
lm(formula = breaks ~ wool * tension, data = warpbreaks)
Multiple Linear Regression 519
Residuals:
Min 1Q Median 3Q Max
-19.5556 -6.8889 -0.6667 7.1944 25.4444
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 44.556 3.647 12.218 2.43e-16 ***
woolB -16.333 5.157 -3.167 0.002677 **
tensionM -20.556 5.157 -3.986 0.000228 ***
tensionH -20.000 5.157 -3.878 0.000320 ***
woolB:tensionM 21.111 7.294 2.895 0.005698 **
woolB:tensionH 10.556 7.294 1.447 0.154327
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 10.94 on 48 degrees of freedom
Multiple R-squared: 0.3778, Adjusted R-squared: 0.3129
F-statistic: 5.828 on 5 and 48 DF, p-value: 0.0002772
Here I've used the cross-factor symbol *, rather than wool + tension +
wool:tension. When both predictors in a two-way interaction are categorical,
there will be a term for each nonreference level of the first predictor combined with all nonreference levels of the second predictor. In this example,
wool is binary with only k = 2 levels and tension has k = 3; therefore, the
only interaction terms present are the "medium" (M) and "high" (H) tension
levels ("low", L, is the reference level) with wool type B (A is the reference
level). Therefore, altogether in the fitted model, there are terms for B, M, H,
B:M, and B:H.
These results provide the same conclusion as the ANOVA analysis-
there is indeed statistical evidence of an interactive effect between wool
type and tension on mean breaks, on top of the contributing main effects
of those predictors.
The general fitted model can be understood as
"Mean warp breaks" = 44.556 ??? 16.333 × "wool type B"
??? 20.556 × "medium tension"
??? 20.000 × "high tension"
+ 21.111 × "wool type B : medium tension"
+ 10.556 × "wool type B : high tension"
The additional interaction terms work the same way as the main
effects-when only categorical predictors are involved, the model can be
seen as a series of additive terms to the overall intercept. Exactly which ones
you use in any given prediction depends on the covariate profile of a given
individual.
520 Chapter 21
Let's have a quick series of examples: for wool A at low tension, the
mean number of warp breaks is predicted as simply the overall intercept; for
wool A at high tension, you have the overall intercept and the main effect
term for high tension; for wool B at low tension, you have the overall intercept and the main effect for wool type B only; and for wool B at medium
tension, you have the overall intercept, the main effect for wool type B, the
main effect for medium tension, and the interactive term for wool B with
medium tension.
You can use predict to estimate the mean warp breaks for these four
scenarios; they're accompanied here with 90 percent confidence intervals:
R> nd <- data.frame(wool=c("A","A","B","B"),tension=c("L","H","L","M"))
R> predict(warp.fit,newdata=nd,interval="confidence",level=0.9)
fit lwr upr
1 44.55556 38.43912 50.67199
2 24.55556 18.43912 30.67199
3 28.22222 22.10579 34.33866
4 28.77778 22.66134 34.89421
21.5.4 Two Continuous
Finally, you'll look at the situation when the two predictors are continuous.
In this case, an interaction term operates as a modifier on the continuous
plane that's fitted using the main effects only. In a similar way to an interaction between a continuous and a categorical predictor, an interaction
between two continuous explanatory variables allows the slope associated
with one variable to be affected, but this time, that modification is made in
a continuous way (that is, according to the value of the other continuous
variable).
Returning to the mtcars data frame, consider MPG once more as a function of horsepower and weight. The fitted model, shown next, includes the
interaction in addition to the main effects of the two continuous predictors.
As you can see, there is a single estimated interactive term, and it is deemed
significantly different from zero.
R> car.fit <- lm(mpg~hp*wt,data=mtcars)
R> summary(car.fit)
Call:
lm(formula = mpg ~ hp * wt, data = mtcars)
Residuals:
Min 1Q Median 3Q Max
-3.0632 -1.6491 -0.7362 1.4211 4.5513
Coefficients:
Estimate Std. Error t value Pr(>|t|)
Multiple Linear Regression 521
(Intercept) 49.80842 3.60516 13.816 5.01e-14 ***
hp -0.12010 0.02470 -4.863 4.04e-05 ***
wt -8.21662 1.26971 -6.471 5.20e-07 ***
hp:wt 0.02785 0.00742 3.753 0.000811 ***
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 2.153 on 28 degrees of freedom
Multiple R-squared: 0.8848, Adjusted R-squared: 0.8724
F-statistic: 71.66 on 3 and 28 DF, p-value: 2.981e-13
The model is written as
"Mean MPG" = 49.80842 ??? 0.12010 × "horsepower"
??? 8.21662 × "weight"
+ 0.02785 × "horsepower : weight"
= 49.80842 ??? 0.12010 × "horsepower"
??? 8.21662 × "weight"
+ 0.02785 × "horsepower" × "weight"
The second version of the model equation provided here reveals for
the first time an interaction expressed as the product of the values of the
two predictors, which is exactly how the fitted model is used to predict the
response. (Technically, this is the same as when at least one of the predictors is categorical-but the dummy coding simply results in zeros and ones
for the respective terms, so multiplication just amounts to the presence or
absence of a given term, as you've seen.)
You can interpret an interaction between two continuous predictors by
considering the sign (+ or ???) of the coefficient. Negativity suggests that as
the values of the predictors increase, the response is reduced after computing the result of the main effects. Positivity, as is the case here, suggests that
as the values of the predictors increase, the effect is an additional increase,
an amplification, on the mean response.
Contextually, the negative main effects of hp and wt indicate that mileage
is naturally reduced for heavier, more powerful cars. However, positivity of
the interactive effect suggests that this impact on the response is "softened"
as horsepower or weight is increased. To put it another way, the negative
relationship imparted by the main effects is rendered "less extreme" as the
values of the predictors get bigger and bigger.
Figure 21-9 contrasts the main-effects-only version of the model
(obtained using lm with the formula mpg~hp+wt; not explicitly fitted in this
section) with the interaction version of the model fitted just above as the
object car.fit.
522 Chapter 21
Figure 21-9: Response surfaces for mean MPG by horsepower and weight, for a maineffects-only model (left), and one that includes the two-way interaction between the
continuous predictors (right)
The plotted response surfaces show the mean MPG on the vertical z-axis
and the two predictor variables on the horizontal axes as marked. You can
interpret the predicted mean MPG, based on a given horsepower and weight
value, as a point on the surface. Note that both surfaces decrease in MPG
(vertically along the z-axis) as you move to larger values of either predictor
along the respective horizontal axes.
I'll show how these plots are created in Chapter 25. For now, they serve
simply to highlight the aforementioned "softening" impact of the interaction in car.fit. On the left, the main-effects-only model shows a flat plane
decreasing according to the negative linear slopes in each predictor. On
the right, however, the presence of the positive interactive term flattens this
plane out, meaning the rate of decrease is slowed as the values of the predictor variables increase.
21.5.5 Higher-Order Interactions
As mentioned, two-way interactions are the most common kind of interactions you'll encounter in applications of regression methods. This is
because for three-way or higher-order terms, you need a lot more data for a
reliable estimation of interactive effects, and there are a number of interpretative complexities to overcome. Three-way interactions are far rarer than
two-way effects, and four-way and above are rarer still.
In Exercise 21.1, you used the nuclear data set found in the boot package (provided with the standard R installation), which includes data on
the constructions of nuclear power plants in the United States. In the exercises, you focused mainly on date and time predictors related to construction permits to model the mean cost of construction for the nuclear power
Multiple Linear Regression 523
plants. For the sake of this example, assume you don't have the data on
these predictors. Can the cost of construction be adequately modeled using
only the variables that describe characteristics of the plant itself?
Load the boot package and access the ?nuclear help page to find details
on the variables: cap (continuous variable describing the capacity of the
plant); cum.n (treated as continuous, describing the number of similar constructions the engineers had previously worked on); ne (binary, describing
whether the plant was in the northeastern United States); and ct (binary,
describing whether the plant had a cooling tower).
The following model is fitted with the final construction cost of the
plant as the response; a main effect for capacity; and main effects of, and all
two-way interactions and the three-way interaction among, cum.n, ne, and ct:
R> nuc.fit <- lm(cost~cap+cum.n*ne*ct,data=nuclear)
R> summary(nuc.fit)
Call:
lm(formula = cost ~ cap + cum.n * ne * ct, data = nuclear)
Residuals:
Min 1Q Median 3Q Max
-162.475 -50.368 -8.833 43.370 213.131
Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept) 138.0336 99.9599 1.381 0.180585
cap 0.5085 0.1127 4.513 0.000157 ***
cum.n -24.2433 6.7874 -3.572 0.001618 **
ne -260.1036 164.7650 -1.579 0.128076
ct -187.4904 76.6316 -2.447 0.022480 *
cum.n:ne 44.0196 12.2880 3.582 0.001577 **
cum.n:ct 35.1687 8.0660 4.360 0.000229 ***
ne:ct 524.1194 200.9567 2.608 0.015721 *
cum.n:ne:ct -64.4444 18.0213 -3.576 0.001601 **
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Residual standard error: 107.3 on 23 degrees of freedom
Multiple R-squared: 0.705, Adjusted R-squared: 0.6024
F-statistic: 6.872 on 8 and 23 DF, p-value: 0.0001264
In this code, you specify the higher-order interactions by extending the
number of variables connected with a * (using * instead of : since you want
to include all the lower-order effects of those three predictors as well).
In the estimated results, the main effect for cap is positive, showing that
an increased power capacity is tied to an increased construction cost. All
other main effects are negative, which at face value seems to imply that a
524 Chapter 21
reduced construction cost is associated with more experienced engineers,
plants constructed in the Northeast, and plants with a cooling tower. However, this isn't an accurate statement since you haven't yet considered the
interactive terms in those predictors. All estimated two-way interactive effects
are positive-having more experienced engineers means a higher construction cost in the Northeast regardless of whether there's a cooling tower, and
having more experienced engineers also means higher costs for plants with a
cooling tower, regardless of region.
Cost is also dramatically increased for plants in the Northeast with a
cooling tower, regardless of the experience of the engineer. All that being
said, the negative three-way interaction suggests that the increased cost associated with more experienced engineers working in the Northeast and on a
plant with a cooling tower is lessened somewhat after the main effects and
two-way interactive effects are calculated.
At the least, this example highlights the complexities associated with
interpreting model coefficients for higher-order interactions. It's also
possible that statistically significant high-order interactions crop up due to
lurking variables that have gone unaccounted for, that is, that the significant
interactions are a spurious manifestation of patterns in the data that simpler
terms involving those missing predictors could explain just as well (if not
better). In part, this motivates the importance of adequate model selection,
which is up next in the discussion.
Exercise 21.3
Return your attention to the cats data frame in package MASS. In the
first few problems in Exercise 21.1, you fitted the main-effect-only
model to predict the heart weights of domestic cats by total body
weight and sex.
a. Fit the model again, and this time include an interaction
between the two predictors. Inspect the model summary.
What do you notice in terms of the parameter estimates and
their significance when compared to the earlier main-effect-only
version?
b. Produce a scatterplot of heart weight on body weight, using
different point characters or colors to distinguish the observations according to sex. Use abline to add two lines denoting
the fitted model. How does this plot differ from the one in
Exercise 21.1 (d)?
c. Predict the heart weight of Tilman's cat using the new model
(remember that Sigma is a 3.4 kg female) accompanied by a
95 percent prediction interval. Compare it to the main-effectsonly model from the earlier exercise.
Multiple Linear Regression 525
In Exercise 21.2, you accessed the trees data frame in the contributed
faraway package. After loading the package, access the ?trees help
file; you'll find the volume and girth measurements you used earlier,
as well as data on the height of each tree.
d. Without using any transformations of the data, fit and inspect a
main-effects-only model for predicting volume from girth and
height. Then, fit and inspect a second version of this model
including an interaction.
e. Repeat (d), but this time use the log transformation of all variables. What do you notice about the significance of the interaction between the untransformed and transformed models? What
does this suggest about the relationships in the data?
Turn back to the mtcars data set and remind yourself of the variables
in the help file ?mtcars.
f. Fit a linear model for mpg based on a two-way interaction between
hp and factor(cyl) and their main effects, as well as a main effect
for wt. Produce a summary of the fit.
g. Interpret the estimated coefficients for the interaction between
horsepower and the (categorical) number of cylinders.
h. Suppose you're keen on purchasing a 1970s performance car.
Your mother advises you to purchase a "practical and economical" car that's capable of an average MPG value of at least 25. You
see three vehicles advertised: car 1 is a four-cylinder, 100 horsepower car that weighs 2100 lbs; car 2 is an eight-cylinder, 210
horsepower car that weighs 3900 lbs; and car 3 is a six-cylinder,
200 horsepower car that weighs 2900 lbs.
i. Use your model to predict the mean MPG for each of the
three cars; provide 95 percent confidence intervals. Based
on your point estimates only, which car would you propose
to your mother?
ii. You still want the most gas-guzzling car you can own with
your mother's blessing, so you decide to be sneaky and
base your decision on what the confidence intervals tell you
instead. Does this change your choice of vehicle?
Important Code in This Chapter
Function/operator Brief description First occurrence
I Include arithmetic term Section 21.4.1, p. 5
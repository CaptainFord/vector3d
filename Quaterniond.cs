using UnityEngine;
using System;
using System.Collections;

namespace UnityEngine {
	public struct Quaterniond {
		public const string defaultNumberFormat = "G5";
		const double dotProductToleranceForEquals = 0.000000001;

		public static void scratchGrounds(){
//			Quaternion q = new Quaternion();
//			UnityTest.QuaterniondTests tests = UnityTest.QuaterniondTests();
			System.Random rand = new System.Random("Nothing at all".GetHashCode());
			Quaterniond qd0 = RandomQuat(rand);
			Quaterniond qd1 = RandomQuat(rand);

			//	Okay, here we go
			Debug.Log ("quats: " + qd0 + " " + qd1);
			for(int i=0; i<4; ++i){
				bool aFirst = i > 1;
				bool inverseA = i % 2 == 0;
				QuaternionDiffTest(qd0, qd1, aFirst, inverseA);
			}

			for(int i=0; i<4; ++i){
				bool aFirst = i > 1;
				bool checkFirst = i % 2 == 0;
				QuaternionDiffTest2(qd0, qd1, aFirst, checkFirst);
			}

		}
		public static void QuaternionDiffTest2(Quaterniond a, Quaterniond b, bool aFirst, bool checkFirst){
			Quaterniond stepOne = aFirst ? a.inverse * b : b * a.inverse;
			string labelOne = aFirst ? "(!a*b)" : "(b*!a)";

			Quaterniond stepTwo = checkFirst ? a * stepOne : stepOne * a;
			string labelTwo = checkFirst ? "a*" + labelOne : labelOne + "*a";
			
			Quaterniond stepThree = checkFirst ? b * stepOne : stepOne * b;
			string labelThree = checkFirst ? "b*" + labelOne : labelOne + "*b";

			Debug.Log (labelTwo + ": " + stepTwo + " - " + b + "\n" + labelThree + ": " + stepThree + " - " + a);
		}
		public static void QuaternionDiffTest(Quaterniond a, Quaterniond b, bool aFirst, bool inverseFirst){
			bool inverseA = inverseFirst ^ !aFirst;
			string labelA = "a", labelB = "b";
			if(inverseA){
				a = Quaterniond.Inverse(a);
				labelA = "!a";
			} else {
				b = Quaterniond.Inverse(b);
				labelB = "!b";
			}
			if(aFirst){
				QuaternionDiffTestFinish(a, b, labelA, labelB);
			} else {
				QuaternionDiffTestFinish(b, a, labelB, labelA);
			}
		}
		public static void QuaternionDiffTestFinish(Quaterniond a, Quaterniond b, string labelA, string labelB){
			Quaterniond result = a * b;
			Debug.Log (labelA + "*" + labelB + "=" + result);
		}
			
			public static Quaterniond RandomQuat(System.Random rand){
			Quaterniond result = new Quaterniond();
			for(int i=0; i<4; ++i){
				result[i] = 2d * rand.NextDouble() - 1d;
			}
			return result;
		}
		public static void scratchGroundsNotRan(){
			Quaternion q = new Quaternion();
			Quaternion q2 = new Quaternion();
			Vector3 v3 = new Vector3();
			Object obj = new Object();
			float f = 0f;

			Debug.Log (Quaternion.identity);
			Debug.Log (Quaternion.identity.w);
			
			q.RegisterInterest(obj);	//	Unnecessary
			q.UnregisterInterest(obj);	//	Unnecessary
			q.GetHashCode();			//	Unnecessary (probably)
			q.SerializeToString();		//	Unnecessary
			q.SerializeXml();			//	Unnecessary

			q.Equals(q2);		//	Unnecessary (but implemented one for Quaternion)
			q.ToString();		//	Implemented
			q.Set(0f,0f,0f,1f);	//	Implemented

			v3 = q.eulerAngles;
			q.SetFromToRotation(v3, v3);
			q.SetLookRotation(v3);
			q.SetLookRotation(v3, v3);
			q.ToAngleAxis(out f, out v3);

//			q = Quaternion.Merge(v3);	//	Not documented. Hmm... wait, that's an extension method from Radical. Ignore it.

			//	Static methods - Mostly constructors
			f = Quaternion.Angle(q, q2);
			q = Quaternion.AngleAxis(f, v3);
			f = Quaternion.Dot(q, q2);
			q = Quaternion.Euler(v3);
			q = Quaternion.FromToRotation(v3, v3);
//			q = Quaternion.identity;	//	Implemented
			q2 = Quaternion.Inverse(q);
			q = Quaternion.Lerp(q, q2, f);
			q = Quaternion.LookRotation(v3);
			q = Quaternion.RotateTowards(q, q2, f);
			q = Quaternion.Slerp(q, q2, f);
		}

		public static Quaterniond identity {
			get {
				return new Quaterniond(0d,0d,0d,1d);
			}
		}

		public static Quaterniond toDoublePrecision(Quaternion q){
			return new Quaterniond(q.x, q.y, q.z, q.w);
		}

		public double x, y, z, w;

		public Quaterniond (double x, double y, double z, double w)
		{
			this.x = x;
			this.y = y;
			this.z = z;
			this.w = w;
		}

		public Quaterniond (Quaternion q)
		{
			this.x = q.x;
			this.y = q.y;
			this.z = q.z;
			this.w = q.w;
		}
		
		public override string ToString(){
			return ToString(defaultNumberFormat);
		}

		public string ToString(string numberFormat){
			return "(" + this.x.ToString(numberFormat) + ", " + 
					this.y.ToString(numberFormat) + ", " + 
					this.z.ToString(numberFormat) + ", " + 
					this.w.ToString(numberFormat) + ")";
		}

		public bool Equals(Quaternion floatQuaternion){
			return this.x == floatQuaternion.x 
					&& this.y == floatQuaternion.y 
					&& this.z == floatQuaternion.z 
					&& this.w == floatQuaternion.w;
		}

		public Vector3d eulerAngles {
			get {
				//	TODO Implement Quaterniond methods
				throw new UnityException("Not Yet Implemented");
			}
		}

		public void Set (double x, double y, double z, double w)
		{
			this.x = x;
			this.y = y;
			this.z = z;
			this.w = w;
		}

		public void SetLookRotation(Vector3d direction){
			SetLookRotation (direction, Vector3d.up);
		}
		public void SetLookRotation(Vector3d direction, Vector3d up){
			//	TODO Implement Quaterniond methods
			throw new UnityException("Not Yet Implemented");
		}

		public void SetFromToRotation(Vector3d from, Vector3d to){
			from = from.normalized;
			to = to .normalized;
			Vector3d cross = Vector3d.Cross(from, to);
			this.x = cross.x;
			this.y = cross.y;
			this.z = cross.z;
			this.w = Mathd.Sqrt(from.sqrMagnitude * to.sqrMagnitude) + Vector3d.Dot(from, to);
			this.Normalize();
		}

		public void ToAngleAxis(out double angle, out Vector3d axis){
			//	http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/
			this.Normalize();

			double halfAngleInRadians = Mathd.Acos(this.w);
			double inverseSinAngle = 1.0d / Mathd.Sqrt (1d - this.w * this.w); 	
//			double inverseSinAngle = 1.0d / Mathd.Sin(halfAngleInRadians);
			angle = halfAngleInRadians * 2.0d * Mathd.Rad2Deg;
			axis = new Vector3d(this.x * inverseSinAngle, 
			                    this.y * inverseSinAngle, 
			                    this.z * inverseSinAngle).normalized;

			//	Huh, these are the same numbers
//			Debug.Log (Mathd.Sqrt(1d - this.w * this.w) + " / " + Mathd.Sin(halfAngleInRadians));

			//	Note: According the webpage above, singularities occur at w=-1 and w=1, 
			//	when the angle is 0 or 180 degrees.
			//	I don't know if this affects the calculation at all.
			//	I appear to be getting non-normalized vectors

//			double angleWouldBe = Mathd.Deg2Rad * angle * 0.5d;
//			Vector3d axisWouldBe = axis.normalized;
//
//			Quaterniond result = new Quaterniond();
//			double sinTheta = Mathd.Sin (angleWouldBe);
//			result.w = Mathd.Cos (angleWouldBe);
//			result.x = axisWouldBe.x * sinTheta;
//			result.y = axisWouldBe.y * sinTheta;
//			result.z = axisWouldBe.z * sinTheta;
//
//			Debug.Log ("angle=" + angle + ";axis=" + axis + ";half= " + halfAngleInRadians + "=>" + angleWouldBe + ";" + axisWouldBe
//			           + "\nquat=" + this + " => " + result + ";sinTheta=" + Mathd.Sin (halfAngleInRadians) + " =>" + sinTheta + ";sanityCheck=" + (inverseSinAngle * sinTheta));
		}

		public Quaterniond inverse {
			get {
				return Quaterniond.Inverse(this);
			}
		}

		public static double Angle(Quaterniond a, Quaterniond b){
			//	The calculation of the angle is part of how ToAngleAxis works
			double angle = Mathd.Acos((a.inverse * b).w) * 2.0d * Mathd.Rad2Deg;
			return angle > 180d ? 360d - angle : angle;
		}

		public static Quaterniond AngleAxis(double angle, Vector3d axis){
			axis = axis.normalized;
			angle = Mathd.Deg2Rad * angle * 0.5d;

			Quaterniond result = new Quaterniond();
			double sinTheta = Mathd.Sin (angle);
			result.w = Mathd.Cos (angle);
			result.x = axis.x * sinTheta;
			result.y = axis.y * sinTheta;
			result.z = axis.z * sinTheta;

			return result;
		}
		public static double Dot(Quaterniond a, Quaterniond b){
			return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
		}
		//	"Returns a rotation that rotates z degrees around the z axis, 
		//	x degrees around the x axis, and y degrees around the y axis (in that order)."
		public static Quaterniond Euler(Vector3d euler){
			return Euler(euler.x, euler.y, euler.z);
		}
		//	"Returns a rotation that rotates z degrees around the z axis, 
		//	x degrees around the x axis, and y degrees around the y axis (in that order)."
		public static Quaterniond Euler(double x, double y, double z){
			return Quaterniond.AngleAxis(y, Vector3d.up)
				* Quaterniond.AngleAxis(x, Vector3d.right) 
					* Quaterniond.AngleAxis(z, Vector3d.forward);
		}

		public static Quaterniond FromToRotation(Vector3d from, Vector3d to){
			Quaterniond result = new Quaterniond();
			result.SetFromToRotation(from, to);
			return result;
		}
		public static Quaterniond FromToRotation(Quaterniond from, Quaterniond to){
			return Quaterniond.Inverse(from) * to;
		}
		public static Quaterniond Inverse(Quaterniond rotation){
			Quaterniond result = new Quaterniond();
			result.x = -rotation.x;
			result.y = -rotation.y;
			result.z = -rotation.z;
			result.w = rotation.w;
			return result;
		}
		public static Quaterniond Lerp(Quaterniond from, Quaterniond to, double t){
			//	TODO Implement Quaterniond methods
			throw new UnityException("Not Yet Implemented");
		}
		public static Quaterniond LookRotation(Vector3d forward){
			return LookRotation (forward, Vector3d.up);
		}
		public static Quaterniond LookRotation(Vector3d forward, Vector3d upwards){
			Quaterniond rotation = new Quaterniond();
			rotation.SetLookRotation(forward, upwards);
			return rotation;
		}
		public static Quaterniond RotateTowards(Quaterniond from, Quaterniond to, double maxDegreesDelta){
			//	TODO Implement Quaterniond methods
			throw new UnityException("Not Yet Implemented");
		}
		public static Quaterniond Slerp(Quaterniond from, Quaterniond to, double t){
			//	TODO Implement Quaterniond methods
			throw new UnityException("Not Yet Implemented");
		}

		public static bool operator !=(Quaterniond lhs, Quaterniond rhs){
			return !(lhs == rhs);
		}
		public static bool operator ==(Quaterniond lhs, Quaterniond rhs){
			//	I think I read something about testing for equivalence basically 
			//	amounting to checking if the dot product was close to zero
			double dotProduct = Quaterniond.Dot (lhs, rhs);
			return Mathd.Abs(dotProduct) < dotProductToleranceForEquals;
		}

		public static Quaterniond operator *(Quaterniond lhs, Quaterniond rhs){
			Quaterniond result = new Quaterniond();
			result.w = lhs.w * rhs.w - lhs.x*rhs.x - lhs.y*rhs.y - lhs.z*rhs.z;
			result.x = lhs.w * rhs.x + lhs.x*rhs.w + lhs.y*rhs.z - lhs.z*rhs.y;
			result.y = lhs.w * rhs.y - lhs.x*rhs.z + lhs.y*rhs.w + lhs.z*rhs.x;
			result.z = lhs.w * rhs.z + lhs.x*rhs.y - lhs.y*rhs.x + lhs.z*rhs.w;
			return result;
		}

		public static Vector3d operator *(Quaterniond lhs, Vector3d rhs){
			return lhs.Multiply(rhs);
		}
		Vector3d Multiply(Vector3d rhs){
			//	Rotation Matrix
			double[,] rotationMatrix = new double[,] {
				{1 - 2*y*y - 2*z*z,	2*x*y - 2*z*w,	2*x*z + 2*y*w},
				{2*x*y + 2*z*w, 1 - 2*x*x - 2*z*z,	2*y*z - 2*x*w},
				{2*x*z - 2*y*w, 2*y*z + 2*x*w,	1 - 2*x*x - 2*y*y}};
			//	3x3 x 3x1 = 3x1	(it's probably this one? Yes)
			//	OR 1x3 x 3x3 = 1x3
			Vector3d result = new Vector3d(
				//	ABij = A
				rhs.x * rotationMatrix[0,0] + rhs.y * rotationMatrix[0,1] + rhs.z * rotationMatrix[0,2],
				rhs.x * rotationMatrix[1,0] + rhs.y * rotationMatrix[1,1] + rhs.z * rotationMatrix[1,2],
				rhs.x * rotationMatrix[2,0] + rhs.y * rotationMatrix[2,1] + rhs.z * rotationMatrix[2,2]
				);
			return result;
		}

		public static Quaterniond operator *(Quaterniond lhs, Quaternion rhs){
			return lhs * new Quaterniond(rhs);
		}
		
		public static Vector3d operator *(Quaterniond lhs, Vector3 rhs){
			return lhs * new Vector3d(rhs);
		}

		public double this[int index]{
			get {
				switch(index){
				case 0:
					return x;
				case 1:
					return y;
				case 2:
					return z;
				case 3:
					return w;
				default:
					throw new ArgumentOutOfRangeException("index: " + index + "  size: 4");
				}
			}
			set {
				switch(index){
				case 0:
					x = value;
					break;
				case 1:
					y = value;
					break;
				case 2:
					z = value;
					break;
				case 3:
					w = value;
					break;
				default:
					throw new ArgumentOutOfRangeException("index: " + index + "  size: 4");
				}
			}


		}

		public void Normalize(){
			double sum = 0;
			for (int i = 0; i < 4; ++i){
				sum += this[i] * this[i];
			}
			double magnitudeInverse = 1d / System.Math.Sqrt(sum);
			for (int i = 0; i < 4; ++i){
				this[i] *= magnitudeInverse;
			}
		}
		public Quaterniond normalized
		{
			get {
				double sum = 0;
				for (int i = 0; i < 4; ++i){
					sum += this[i] * this[i];
				}
				double magnitudeInverse = 1d / System.Math.Sqrt(sum);
				return new Quaterniond(this.x * magnitudeInverse,
				                       this.y * magnitudeInverse,
				                       this.z * magnitudeInverse,
				                       this.w * magnitudeInverse);
			}
		}

		public double GetSumOfSquares(){
			double sum = 0;
			for(int i=0; i<4; ++i){
				sum += this[i] * this[i];
			}
			return sum;
		}

	}
}

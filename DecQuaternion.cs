using UnityEngine;
using System;
using System.Collections;

namespace UnityEngine {
	public struct DecQuaternion {
		public const string defaultNumberFormat = "G5";

		public static DecQuaternion identity {
			get {
				return new DecQuaternion(0m,0m,0m,1m);
			}
		}

		public decimal x, y, z, w;

		public DecQuaternion (decimal x, decimal y, decimal z, decimal w)
		{
			this.x = x;
			this.y = y;
			this.z = z;
			this.w = w;
		}

		public DecQuaternion (double x, double y, double z, double w)
		{
			this.x = (decimal)x;
			this.y = (decimal)y;
			this.z = (decimal)z;
			this.w = (decimal)w;
		}
		public DecQuaternion (float x, float y, float z, float w) 
				: this((decimal)x,(decimal)y,(decimal)z,(decimal)w){}
		public DecQuaternion (Quaterniond q) 
				: this((decimal)q.x, (decimal)q.y, (decimal)q.z, (decimal)q.w) {}
		public DecQuaternion (Quaternion q) 
				: this((decimal)q.x, (decimal)q.y, (decimal)q.z, (decimal)q.w) {}
		
		public override string ToString(){
			return ToString(defaultNumberFormat);
		}

		public string ToString(string numberFormat){
			return "(" + this.x.ToString(numberFormat) + ", " + 
					this.y.ToString(numberFormat) + ", " + 
					this.z.ToString(numberFormat) + ", " + 
					this.w.ToString(numberFormat) + ")";
		}

		public override bool Equals(object obj){
			if(!(obj is DecQuaternion)){
				return false;
			}
			DecQuaternion other = (DecQuaternion)obj;
			return this.x.Equals(other.x)
					&& this.y.Equals(other.y)
					&& this.z.Equals(other.z)
					&& this.w.Equals(other.w);
		}

		public override int GetHashCode(){
			return this.x.GetHashCode() 
					^ this.y.GetHashCode() << 2 
					^ this.z.GetHashCode() << 4 
					^ this.w.GetHashCode() >> 2;
		}
		public bool EqualTo(Quaterniond doubleQuaternion){
			return this.x.Equals(doubleQuaternion.x)
				&& this.y.Equals(doubleQuaternion.y)
					&& this.z.Equals(doubleQuaternion.z)
					&& this.w.Equals(doubleQuaternion.w);
		}
		public bool EqualTo(Quaternion floatQuaternion){
			return this.x.Equals(floatQuaternion.x)
					&& this.y.Equals(floatQuaternion.y)
	                && this.z.Equals(floatQuaternion.z)
	                && this.w.Equals(floatQuaternion.w);
		}

		public void Set (decimal x, decimal y, decimal z, decimal w)
		{
			this.x = x;
			this.y = y;
			this.z = z;
			this.w = w;
		}

		public void SetLookRotation(DecVector3 direction){
			//	Okay, what we actually have going on here is two rotations. A heading and a pitch.
			//	The heading is a rotation around the given up axis. The pitch, then, must be around an axis perpendicular to the up direction and the heading.
			//	Not perpendicular to the direction vector and up, as I've seen in other algorithms, but perpendicular to the heading and the up axis.
			DecVector3 up = DecVector3.up;
			DecVector3 right = DecVector3.Cross(up, direction);    // The perpendicular vector to Up and Direction
			DecVector3 heading = DecVector3.Cross(right, up);   
			this = DecQuaternion.FromToRotation(heading, direction) * DecQuaternion.FromToRotation(DecVector3.forward, heading);
//			SetLookRotation(direction, DecVector3.up);
		}
		public void SetLookRotation(DecVector3 direction, DecVector3 up){
			//	Okay, after some experimenting, I'm noticing that this version is rotating twice as far as it should when varying the up vector.
			//	Such that when the up vector points straight down, it's the same as pointing straight up.

			//	Ultimately, I think the problem stems from the fact that I assumed the camera would start tilted, which it doesn't

			DecVector3 right = DecVector3.Cross(up, direction);    // The perpendicular vector to Up and Direction
			DecVector3 heading = DecVector3.Cross(right, up);   
			this = DecQuaternion.FromToRotation(heading, direction) * DecQuaternion.FromToRotation(DecVector3.forward, heading);

			//	This correctly orients the camera along the desired axis. Now it just needs to pitch from there to ... well ... that.
			this = DecQuaternion.FromToRotation(DecVector3.up, up);

			this = DecQuaternion.FromToRotation(heading, direction) * DecQuaternion.FromToRotation(this * DecVector3.forward, heading) * this;
			// And ... that's a wrap.
		}

		public void SetFromToRotation(DecVector3 from, DecVector3 to){
			from = from.normalized;
			to = to .normalized;
			DecVector3 axis = DecVector3.Cross(from, to);
			this.x = axis.x;
			this.y = axis.y;
			this.z = axis.z;
			this.w = MathDec.Sqrt(from.sqrMagnitude * to.sqrMagnitude) + DecVector3.Dot(from, to);
			this.Normalize();
		}

		public void ToAngleAxis(out decimal angle, out DecVector3 axis){
			//	http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/
			this.Normalize();

			decimal halfAngleInRadians = MathDec.Acos(this.w);
			decimal inverseSinAngle = 1.0m / MathDec.Sqrt (1m - this.w * this.w);
			angle = halfAngleInRadians * 2.0m * MathDec.Rad2Deg;
			axis = new DecVector3(this.x * inverseSinAngle, 
			                    this.y * inverseSinAngle, 
			                    this.z * inverseSinAngle).normalized;
		}

		public DecQuaternion inverse {
			get {
				return DecQuaternion.Inverse(this);
			}
		}

		public static decimal Angle(DecQuaternion a, DecQuaternion b){
			//	The calculation of the angle is part of how ToAngleAxis works
			decimal angle = MathDec.Acos((a.inverse * b).w) * 2.0m * MathDec.Rad2Deg;
			return angle > 180m ? 360m - angle : angle;
		}

		public static DecQuaternion AngleAxis(decimal angle, DecVector3 axis){
			axis = axis.normalized;
			angle = MathDec.Deg2Rad * angle * 0.5m;

			DecQuaternion result = new DecQuaternion();
			decimal sinTheta = MathDec.Sin (angle);
			result.w = MathDec.Cos (angle);
			result.x = axis.x * sinTheta;
			result.y = axis.y * sinTheta;
			result.z = axis.z * sinTheta;

			return result;
		}
		public static decimal Dot(DecQuaternion a, DecQuaternion b){
			return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
		}

		public static DecQuaternion FromToRotation(DecVector3 from, DecVector3 to){
			DecQuaternion result = new DecQuaternion();
			result.SetFromToRotation(from, to);
			return result;
		}
		public static DecQuaternion FromToRotation(DecQuaternion from, DecQuaternion to){
			return DecQuaternion.Inverse(from) * to;
		}
		public static DecQuaternion Inverse(DecQuaternion rotation){
			DecQuaternion result = new DecQuaternion();
			result.x = -rotation.x;
			result.y = -rotation.y;
			result.z = -rotation.z;
			result.w = rotation.w;
			return result;
		}
		public static DecQuaternion Lerp(DecQuaternion from, DecQuaternion to, decimal t){
			from = from.normalized;
			to = to.normalized;
			decimal multA = 1 - t;
			decimal multB = t;

			if(DecQuaternion.Dot(from, to) < 0m){
				from = -from;
			}
			
			return new DecQuaternion(
				from.x * multA + to.x * multB,
				from.y * multA + to.y * multB,
				from.z * multA + to.z * multB,
				from.w * multA + to.w * multB
				).normalized;
		}
		public static DecQuaternion LookRotation(DecVector3 forward){
			return LookRotation (forward, DecVector3.up);
		}
		public static DecQuaternion LookRotation(DecVector3 forward, DecVector3 upwards){
			DecQuaternion rotation = new DecQuaternion();
			rotation.SetLookRotation(forward, upwards);
			return rotation;
		}
		public static DecQuaternion RotateTowards(DecQuaternion from, DecQuaternion to, decimal maxDegreesDelta){
			decimal angle = DecQuaternion.Angle(from, to);
			decimal t = angle > maxDegreesDelta  ? maxDegreesDelta / angle : 1m;
//			DecQuaternion difference = to.inverse * from;
			return DecQuaternion.Slerp(from, to, t);
		}
		public static DecQuaternion Slerp(DecQuaternion from, DecQuaternion to, decimal t){
			from = from.normalized;
			to = to.normalized;

			decimal dot = DecQuaternion.Dot(from, to);
			if(dot < 0m){
				from = -from;
				dot = -dot;
			}
			if(dot > 0.9995m){
				return Lerp (from, to, t);
			}
//			dot = MathDec.Clamp(dot, -1m, 1m);

			decimal theta =  MathDec.Acos(dot) * t;
			;
			decimal multA = MathDec.Cos(theta);
			decimal multB = MathDec.Sin(theta);

			DecQuaternion c = to - from * dot;
			c.Normalize();

			return new DecQuaternion(
				from.x * multA + c.x * multB,
				from.y * multA + c.y * multB,
				from.z * multA + c.z * multB,
				from.w * multA + c.w * multB
				).normalized;
		}

		public static DecQuaternion operator -(DecQuaternion q){
			return new DecQuaternion(-q.x, -q.y, -q.z, -q.w);
		}

		public static bool operator !=(DecQuaternion lhs, DecQuaternion rhs){
			return !(lhs == rhs);
		}
		public static bool operator ==(DecQuaternion lhs, DecQuaternion rhs){
			return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w;
		}

		public static DecQuaternion operator *(DecQuaternion lhs, DecQuaternion rhs){
			DecQuaternion result = new DecQuaternion();
			result.w = lhs.w * rhs.w - lhs.x*rhs.x - lhs.y*rhs.y - lhs.z*rhs.z;
			result.x = lhs.w * rhs.x + lhs.x*rhs.w + lhs.y*rhs.z - lhs.z*rhs.y;
			result.y = lhs.w * rhs.y - lhs.x*rhs.z + lhs.y*rhs.w + lhs.z*rhs.x;
			result.z = lhs.w * rhs.z + lhs.x*rhs.y - lhs.y*rhs.x + lhs.z*rhs.w;
			return result;
		}

		public static DecQuaternion operator *(DecQuaternion lhs, decimal rhs){
			DecQuaternion result = new DecQuaternion();
			result.w = lhs.w * rhs;
			result.x = lhs.x * rhs;
			result.y = lhs.y * rhs;
			result.z = lhs.z * rhs;
			return result;
		}
		
		public static DecQuaternion operator *(DecQuaternion lhs, Quaternion rhs){
			return lhs * new DecQuaternion(rhs);
		}

		public static DecQuaternion operator +(DecQuaternion lhs, DecQuaternion rhs){
			return new DecQuaternion(lhs.x + rhs.x,
			                       lhs.y + rhs.y, 
			                       lhs.z + rhs.z, 
			                       lhs.w + rhs.w);
		}	
		public static DecQuaternion operator -(DecQuaternion lhs, DecQuaternion rhs){
			return new DecQuaternion(lhs.x - rhs.x,
			                       lhs.y - rhs.y, 
			                       lhs.z - rhs.z, 
			                       lhs.w - rhs.w);
		}

		public static DecVector3 operator *(DecQuaternion lhs, Vector3 rhs){
			return lhs.Multiply(rhs);
		}
		public static DecVector3 operator *(DecQuaternion lhs, DecVector3 rhs){
			return lhs.Multiply(rhs);
		}
		public DecVector3 Multiply(Vector3 rhs){
			return Multiply ((DecVector3)rhs);
		}
		public DecVector3 Multiply(DecVector3 rhs){
			decimal xx = x*x*2;
			decimal xy = x*y*2;
			decimal xz = x*z*2;
			decimal xw = x*z*2;
			decimal yy = y*y*2;
			decimal yz = y*z*2;
			decimal yw = y*z*2;
			decimal zz = z*z*2;
			decimal zw = z*z*2;
			return new DecVector3(
				//	ABij = A
				rhs.x * (1 - yy - zz) 	+ rhs.y * (xy - zw) 	+ rhs.z * (xz + yw),
				rhs.x * (xy + zw) 		+ rhs.y * (1 - xx - zz) + rhs.z * (yz - xw),
				rhs.x * (xz - yw) 		+ rhs.y * (yz + xw) 	+ rhs.z * (1 - xx - yy)
				);
		}
		public DecVector3 InverseMultiply(Vector3 rhs){
			return InverseMultiply((DecVector3)rhs);
		}
		public DecVector3 InverseMultiply(DecVector3 rhs){
			decimal xx = x*x*2;
			decimal xy = x*y*2;
			decimal xz = x*z*2;
			decimal xw = -x*z*2;
			decimal yy = y*y*2;
			decimal yz = y*z*2;
			decimal yw = -y*z*2;
			decimal zz = z*z*2;
			decimal zw = -z*z*2;
			return new DecVector3(
				//	ABij = A
				rhs.x * (1 - yy - zz) 	+ rhs.y * (xy - zw) 	+ rhs.z * (xz + yw),
				rhs.x * (xy + zw) 		+ rhs.y * (1 - xx - zz) + rhs.z * (yz - xw),
				rhs.x * (xz - yw) 		+ rhs.y * (yz + xw) 	+ rhs.z * (1 - xx - yy)
				);
		}
//		DecVector3 MultiplyOld(DecVector3 rhs){
//			decimal[,] rotationMatrix = new decimal[,] {
//				{1 - 2*y*y - 2*z*z,	2*x*y - 2*z*w,	2*x*z + 2*y*w},
//				{2*x*y + 2*z*w, 1 - 2*x*x - 2*z*z,	2*y*z - 2*x*w},
//				{2*x*z - 2*y*w, 2*y*z + 2*x*w,	1 - 2*x*x - 2*y*y}};
//			DecVector3 result = new DecVector3(
//				//	ABij = A
//				rhs.x * rotationMatrix[0,0] + rhs.y * rotationMatrix[0,1] + rhs.z * rotationMatrix[0,2],
//				rhs.x * rotationMatrix[1,0] + rhs.y * rotationMatrix[1,1] + rhs.z * rotationMatrix[1,2],
//				rhs.x * rotationMatrix[2,0] + rhs.y * rotationMatrix[2,1] + rhs.z * rotationMatrix[2,2]
//				);
//			return result;
//		}


		public decimal this[int index]{
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
			decimal sum = this.w * this.w + this.x * this.x + this.y * this.y + this.z * this.z;
//			decimal magnitudeInverse = MathDec.Sign (this.w) / System.Math.Sqrt(sum);
			decimal magnitudeInverse = 1m / MathDec.Sqrt(sum);
			this.w *= magnitudeInverse;
			this.x *= magnitudeInverse;
			this.y *= magnitudeInverse;
			this.z *= magnitudeInverse;
		}
//		public void NormalizeOld(){
//			decimal sum = 0;
//			for (int i = 0; i < 4; ++i){
//				sum += this[i] * this[i];
//			}
//			decimal magnitudeInverse = MathDec.Sign (this.w) / System.Math.Sqrt(sum);
//			for (int i = 0; i < 4; ++i){
//				this[i] *= magnitudeInverse;
//			}
//		}
		public void NormalizeSign(){
			if(this.w < 0m){
				this.w = -this.w;
				this.x = -this.x;
				this.y = -this.y;
				this.z = -this.z;
			}
		}

		public DecQuaternion normalized
		{
			get {
				decimal sum = 0;
				for (int i = 0; i < 4; ++i){
					sum += this[i] * this[i];
				}
				decimal magnitudeInverse = 1m / MathDec.Sqrt(sum);
				return new DecQuaternion(this.x * magnitudeInverse,
				                       this.y * magnitudeInverse,
				                       this.z * magnitudeInverse,
				                       this.w * magnitudeInverse);
			}
		}

		public decimal GetSumOfSquares(){
			decimal sum = 0;
			for(int i=0; i<4; ++i){
				sum += this[i] * this[i];
			}
			return sum;
		}

		public Quaternion ToFloatQuaternion(){
			return new Quaternion((float)x,(float)y,(float)z,(float)w);
		}

		/**
		 *	Tests if the values of the vector are approximately equal. The value of precisionBase is a multiple of the
		 *	base precision of a decimal, which is itself multiplied by the values to determine the tolerance. Any value
		 *	less than zero means "must be exactly zero".
		 */
		public static bool Approximately(DecQuaternion lhs, DecQuaternion rhs, decimal toleranceMultiplier){
			return MathDec.ApproximatelyRelative(lhs.x, rhs.x, toleranceMultiplier) 
				&& MathDec.ApproximatelyRelative(lhs.y, rhs.y, toleranceMultiplier)
					&& MathDec.ApproximatelyRelative(lhs.z, rhs.z, toleranceMultiplier)
			&& MathDec.ApproximatelyRelative(lhs.w, rhs.w, toleranceMultiplier);
		}
		
		public static bool ApproximatelyFixed(DecQuaternion lhs, DecQuaternion rhs, decimal toleranceBase){
			return MathDec.ApproximatelyFixed(lhs.x, rhs.x, toleranceBase) 
					&& MathDec.ApproximatelyFixed(lhs.y, rhs.y, toleranceBase)
					&& MathDec.ApproximatelyFixed(lhs.z, rhs.z, toleranceBase)
					&& MathDec.ApproximatelyFixed(lhs.w, rhs.w, toleranceBase);
		}

		static public implicit operator DecQuaternion(Quaternion value) {
			return new DecQuaternion(value);
		}
		static public explicit operator Quaternion(DecQuaternion value)
		{
			return value.ToFloatQuaternion();
		}
	}
}

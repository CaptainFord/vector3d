using UnityEngine;
using System;
using System.Collections;

namespace UnityEngine {
	public struct Quaterniond {
		public const string defaultNumberFormat = "G5";

		public static void scratchGrounds(){
//			Quaternion q = new Quaternion();
//			UnityTest.QuaterniondTests tests = UnityTest.QuaterniondTests();
//			System.Random rand = new System.Random("Nothing at all".GetHashCode());
//			Quaterniond qd0 = RandomQuat(rand);
//			Quaterniond qd1 = RandomQuat(rand);

//			Debug.Log ("quats: " + qd0 + " " + qd1);
//			for(int i=0; i<4; ++i){
//				bool aFirst = i > 1;
//				bool inverseA = i % 2 == 0;
//				QuaternionDiffTest(qd0, qd1, aFirst, inverseA);
//			}
//
//			for(int i=0; i<4; ++i){
//				bool aFirst = i > 1;
//				bool checkFirst = i % 2 == 0;
//				QuaternionDiffTest2(qd0, qd1, aFirst, checkFirst);
//			}

//			Debug.Log (Convert.ToString(BitConverter.DoubleToInt64Bits(0.0),2) + "\n" 
//			           + Convert.ToString(BitConverter.DoubleToInt64Bits(0.0 / 1.0),2));
//			Vector3dEqualsInvestigation();
//			QuaternionEqualsInvestigation();

			Quaternion q = new Quaternion(-0.4850f, 0.3952f, -0.6186f, -0.4753f);
			Debug.Log ("q.eulerAngles=" + q.eulerAngles + "\n" + new Quaterniond(q).eulerAngles);
		}
		static void Vector3dEqualsInvestigation(){
			TestVector3Equals(0,0,0);
			TestVector3Equals(1,1,1);

		}
		static void TestVector3Equals(double x, double y, double z){
			Vector3d dv0 = new Vector3d(x,y,z);
			Vector3d dv1 = new Vector3d(x,y,z);
			Vector3 fv0 = (Vector3)dv0;
			Vector3 fv1 = (Vector3)dv1;

			Debug.Log (EqualsMessage(dv0, dv1) + "\n" + EqualsMessage(fv0, fv1));


		}
		static string EqualsMessage(Vector3d a, Vector3d b){
			return (a == b) + " = (" + a + " == " + b + ") = " + (a.Equals(b));
		}
		static string EqualsMessage(Vector3 a, Vector3 b){
			return (a == b) + " = (" + a + " == " + b + ") = " + (a.Equals(b));
		}
		static void LogEquals(Vector3 a, Vector3 b){
			Debug.Log ((a == b) + " = (" + a + " == " + b + ")");
		}
		static void LogEquals(object a, object b){
			Debug.Log ((a == b) + " = (" + a + " == " + b + ")");
		}
		
		static void QuaternionEqualsInvestigation(){
//			QuaternionEqualsInvestigation(1f, 1f, 1f, 1f);
//			Quaternion minEquivalent;
//			Quaternion maxEquivalent;


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

		public override bool Equals(object obj){
			if(!(obj is Quaterniond)){
				return false;
			}
			Quaterniond other = (Quaterniond)obj;
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

		public bool EqualTo(Quaternion floatQuaternion){
			return this.x.Equals(floatQuaternion.x)
					&& this.y.Equals(floatQuaternion.y)
	                && this.z.Equals(floatQuaternion.z)
	                && this.w.Equals(floatQuaternion.w);
		}

		public void Set (double x, double y, double z, double w)
		{
			this.x = x;
			this.y = y;
			this.z = z;
			this.w = w;
		}

		public void SetLookRotation(Vector3d direction){
			//	Okay, what we actually have going on here is two rotations. A heading and a pitch.
			//	The heading is a rotation around the given up axis. The pitch, then, must be around an axis perpendicular to the up direction and the heading.
			//	Not perpendicular to the direction vector and up, as I've seen in other algorithms, but perpendicular to the heading and the up axis.
			Vector3d up = Vector3d.up;
			Vector3d right = Vector3d.Cross(up, direction);    // The perpendicular vector to Up and Direction
			Vector3d heading = Vector3d.Cross(right, up);   
			this = Quaterniond.FromToRotation(heading, direction) * Quaterniond.FromToRotation(Vector3d.forward, heading);
//			SetLookRotation(direction, Vector3d.up);
		}
		public void SetLookRotation(Vector3d direction, Vector3d up){
			//	Okay, after some experimenting, I'm noticing that this version is rotating twice as far as it should when varying the up vector.
			//	Such that when the up vector points straight down, it's the same as pointing straight up.

			//	Ultimately, I think the problem stems from the fact that I assumed the camera would start tilted, which it doesn't

			Vector3d right = Vector3d.Cross(up, direction);    // The perpendicular vector to Up and Direction
			Vector3d heading = Vector3d.Cross(right, up);   
			this = Quaterniond.FromToRotation(heading, direction) * Quaterniond.FromToRotation(Vector3d.forward, heading);

			//	This correctly orients the camera along the desired axis. Now it just needs to pitch from there to ... well ... that.
			this = Quaterniond.FromToRotation(Vector3d.up, up);

			this = Quaterniond.FromToRotation(heading, direction) * Quaterniond.FromToRotation(this * Vector3d.forward, heading) * this;
			// And ... that's a wrap.
		}

		public void SetFromToRotation(Vector3d from, Vector3d to){
			from = from.normalized;
			to = to .normalized;
			Vector3d axis = Vector3d.Cross(from, to);
			this.x = axis.x;
			this.y = axis.y;
			this.z = axis.z;
			this.w = Mathd.Sqrt(from.sqrMagnitude * to.sqrMagnitude) + Vector3d.Dot(from, to);
			this.Normalize();
		}

		public void ToAngleAxis(out double angle, out Vector3d axis){
			//	http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/
			this.Normalize();

			double halfAngleInRadians = Mathd.Acos(this.w);
			double inverseSinAngle = 1.0d / Mathd.Sqrt (1d - this.w * this.w);
			angle = halfAngleInRadians * 2.0d * Mathd.Rad2Deg;
			axis = new Vector3d(this.x * inverseSinAngle, 
			                    this.y * inverseSinAngle, 
			                    this.z * inverseSinAngle).normalized;
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
//		public static Quaterniond EulerOld(double x, double y, double z){
//			return Quaterniond.AngleAxis(y, Vector3d.up)
//					* Quaterniond.AngleAxis(x, Vector3d.right) 
//					* Quaterniond.AngleAxis(z, Vector3d.forward);
//		}
		//	I would like to do some benchmarking to compare the old version to this new one
		//	"Returns a rotation that rotates z degrees around the z axis, 
		//	x degrees around the x axis, and y degrees around the y axis (in that order)."
		public static Quaterniond Euler(double xIn, double yIn, double zIn){
			double mult = Mathd.Deg2Rad * 0.5d;
			double x = xIn * mult;
			double y = yIn * mult;
			double z = zIn * mult;
			
			double aroundX_x, aroundX_w;
			double aroundY_y, aroundY_w;
			double aroundZ_z, aroundZ_w;
			
			{
				aroundY_w = Mathd.Cos (y);
				aroundY_y = Mathd.Sin (y);
			}
			{
				aroundX_w = Mathd.Cos (x);
				aroundX_x = Mathd.Sin (x);
			}
			{
				aroundZ_w = Mathd.Cos (z);
				aroundZ_z = Mathd.Sin (z);
			}
			
			double xy_x, xy_y, xy_z, xy_w;
			double xyz_x, xyz_y, xyz_z, xyz_w;
			
			{
				xy_w =  aroundY_w * aroundX_w;
				xy_x =  aroundY_w * aroundX_x;
				xy_y =  aroundY_y * aroundX_w;
				xy_z = -aroundY_y * aroundX_x;
			}
			{
				xyz_w =  xy_w * aroundZ_w - xy_z * aroundZ_z;
				xyz_x =  xy_x * aroundZ_w + xy_y * aroundZ_z;
				xyz_y = -xy_x * aroundZ_z + xy_y * aroundZ_w;
				xyz_z =  xy_w * aroundZ_z + xy_z * aroundZ_w;
			}
			return new Quaterniond(xyz_x, xyz_y, xyz_z, xyz_w);
		}
		public Vector3d eulerAnglesBlegh{
			get {
				//	ZXYs = 3,1,1,0
				//	ZXYr = 2,0,1,1
				//				int[] tuple = {3,1,1,0};
				return EulerAngles (false, 3,1,1,0);
//				return EulerAngles (2,0,1,1);
//				return EulerAngles (1,1,0,0);
			}
		}
		public Vector3d EulerAngles(bool flipMatrix, params int[] tuple){
//			〈a, b〉 := 〈My,x, Mz,x〉
//			〈c, s, r〉 := Givens(a, b)
//			sx := c*Mz,y - s*My,y
//			cx := c*Mz,z - s*My,z
//			θx := atan2(sx, cx)
//			θz := atan2(r, Mx,x)
//			θx′ := atan2(s, c)
//			i 	n 	a 	r 		  	Sys.   		i 	n 	a 	r 		  	Sys.   		i 	n 	a 	r 		  	Sys.   		i 	n 	a 	r 		  	Sys.
//		〈 	1 	0 	0 	0 	〉 		xzxs 	〈 	1 	0 	1 	0 	〉 		xzys 	〈 	1 	1 	0 	0 	〉 		xyxs 	〈 	1 	1 	1 	0 	〉 		xyzs
//		〈 	2 	0 	0 	0 	〉 		yxys 	〈 	2 	0 	1 	0 	〉 		yxzs 	〈 	2 	1 	0 	0 	〉 		yzys 	〈 	2 	1 	1 	0 	〉 		yzxs
//		〈 	3 	0 	0 	0 	〉 		zyzs 	〈 	3 	0 	1 	0 	〉 		zyxs 	〈 	3 	1 	0 	0 	〉 		zxzs 	〈 	3 	1 	1 	0 	〉 		zxys
//		〈 	1 	0 	0 	1 	〉 		xzxr 	〈 	1 	0 	1 	1 	〉 		yzxr 	〈 	1 	1 	0 	1 	〉 		xyxr 	〈 	1 	1 	1 	1 	〉 		zyxr
//		〈 	2 	0 	0 	1 	〉 		yxyr 	〈 	2 	0 	1 	1 	〉 		zxyr 	〈 	2 	1 	0 	1 	〉 		yzyr 	〈 	2 	1 	1 	1 	〉 		xzyr
//		〈 	3 	0 	0 	1 	〉 		zyzr 	〈 	3 	0 	1 	1 	〉 		xyzr 	〈 	3 	1 	0 	1 	〉 		zxzr 	〈 	3 	1 	1 	1 	〉 		yxzr

			//http://web.archive.org/web/20110722193141/http://cgafaq.info/wiki/Euler_angles_from_matrix
			double[,] m = new double[,] {
				{1 - 2*y*y - 2*z*z,	2*x*y - 2*z*w,	2*x*z + 2*y*w},
				{2*x*y + 2*z*w, 1 - 2*x*x - 2*z*z,	2*y*z - 2*x*w},
				{2*x*z - 2*y*w, 2*y*z + 2*x*w,	1 - 2*x*x - 2*y*y}};

			if(flipMatrix){
				//	flips the matrix
//				m = new double[,] {
//					{m[0,0],m[0,1],m[0,2]},
//					{m[1,0],m[1,1],m[1,2]},
//					{m[2,0],m[2,1],m[2,2]},
//
//				};
				m = new double[,] {
					{m[0,0],m[1,0],m[2,0]},
					{m[0,1],m[1,1],m[2,1]},
					{m[0,2],m[1,2],m[2,2]},
				};
			}

//			int[] EulNext = {y,z,x,y}
			int[] eulNext = {1,2,3,1,2,3};
			int i = tuple[0];
			int neg = tuple[1];
			int alt = tuple[2];
			int rev = tuple[3];

			int j = eulNext[i + neg];
			int k = 6 - i - j;
			int h = eulNext[k+(1^neg^alt)];
			j--;k--;h--;i--;

//			Debug.Log ("i=" + i + " j=" + j + " k=" + k + "h=" + h);
			Vector3d v = new Vector3d(m[0,i],m[1,i],m[2,i]);
			double a = v[h], b = v[k];
			double c, s, r;
			EulerGivens(out c, out s, out r, a, b);
			v[h] = r;

			double s1 = c * m[k,j] - s * m[h,j];
			double c1 = c * m[k,k] - s * m[h,k];

			Vector3d result = new Vector3d();
			result.x = Mathd.Atan2(s1, c1);
			result.y = Mathd.Atan2(v[j], v[i]);
			result.z = Mathd.Atan2(s, c);

			if(alt == 1){
				result.z = -result.z;
			}
			if(neg == 1){
				result = -result;
			}
			if(rev == 1){
				result = new Vector3d(result.z, result.y, result.x);
			}
			for(int l=0; l<3; ++l){
				result[l] *= Mathd.Rad2Deg;
				if(result[l] < 0d){
					result[l] += 360d;
				}
			}
			//	Funny thing, actually, the rotations aren't actually input in the order they're done. 
			//	Instead, they're input by the axis they rotate around.
			//	Since it is defined as z, x, y, and the return is in order performed, they need to be reordered like so:
			result = new Vector3d(result.y, result.z, result.x);
//			result = new Vector3d(result.z, result.y, result.x);
			return result;

//		〈j, k, h〉 := indices(i, neg, alt)
//		v∗ := M∗,i
//		〈a, b〉 := 〈vh, vk〉
//		〈c, s, vh〉 := Givens(a, b)
//		s1 := c*Mk,j - s*Mh,j
//		c1 := c*Mk,k - s*Mh,k
//		θ1 := atan2(s1, c1)
//		θ2 := atan2(vj, vi)
//		θ3 := atan2(s, c)
//		if (alt = 1) then θ3 := -θ3
//		if (neg = 1) then 〈θ1, θ2, θ3〉 := 〈-θ1, -θ2, -θ3〉
//		if (rev = 1) then 〈θ1, θ2, θ3〉 := 〈θ3, θ2, θ1〉
		


		}

		void EulerGivens(out double c, out double s, out double r, double a, double b){
			if(b == 0){
				c = Mathd.Sign(a);
				s = 0;
				r = Mathd.Abs(a);
			} else if(a == 0){
				c = 0;
				s = Mathd.Sign(b);
				r = Mathd.Abs(b);
			} else if(Mathd.Abs(b) > Mathd.Abs(a)){
				double t = ((double)a) / b;
				double u = Mathd.Sign(b)*Mathd.Sqrt(1 + t * t);
				s = 1 / u;
				c = s * t;
				r = b * u;
			} else {
				double t = ((double)b) / a;
				double u = Mathd.Sign(a)*Mathd.Sqrt(1 + t * t);
				c = 1 / u;
				s = c * t;
				r = a * u;
			}
		}
		//	"Returns a rotation that rotates z degrees around the z axis, 
		//	x degrees around the x axis, and y degrees around the y axis (in that order)."
		public Vector3d eulerAnglesZYX {
			get {
				double sqx = x * x;
				double sqy = y * y;
				double sqz = z * z;
				double sqw = w * w;
				
				double unit = sqx + sqy + sqz + sqw; // if normalised is one, otherwise is correction factor
				double test = x * y + w * z;
				Vector3d v;
				
				if (test>0.49999995d*unit) { // singularity at north pole
					v.y = 2d * Mathd.Atan2 (y, x) * Mathd.Rad2Deg;
					v.x = 90d;
					v.z = 0d;
					return v;
				}
				if (test<-0.49999995d*unit) { // singularity at south pole
					v.y = -2d * Mathd.Atan2 (y, x) * Mathd.Rad2Deg;
					v.x = -90d;
					v.z = 0d;
					return v;
				}
				v.y = Mathd.Rad2Deg * Mathd.Atan2 (2 * (z*y + w*x), 1 - 2 * (sqx + sqy));     // Yaw
				v.x = Mathd.Rad2Deg * Mathd.Asin (-2 * (x*z + w*y));                             // Pitch
				v.z = Mathd.Rad2Deg * Mathd.Atan2 (2 * (w*z + x*y), 1 - 2 * (sqz + sqx));      // Roll
				
				for(int i=0; i<3; ++i){
					if(v[i] < 0d){
						v[i] += 360d;
					}
				}
				return v;

			}
		}
		//	"Returns a rotation that rotates z degrees around the z axis, 
		//	x degrees around the x axis, and y degrees around the y axis (in that order)."
		public Vector3d eulerAngles {
			get {
				//	I learned a lot from optimizing the Euler method
				//	The most critical thing I learned is that the three quaternions generated
				//	by Euler all have only two components. Each vector component belongs to only
				//	one of the quaternions, and all three have a real component.

				double sqx = x * x;
				double sqy = y * y;
				double sqz = z * z;
				double sqw = w * w;

				double unit = sqx + sqy + sqz + sqw; // if normalised is one, otherwise is correction factor
				double test = x * w - y * z;
				Vector3d v;
				
				if (test>0.49999995d*unit) { // singularity at north pole
					v.y = 2d * Mathd.Atan2 (y, x) * Mathd.Rad2Deg;
					v.x = 90d;
					v.z = 0d;
					return v;
				}
				if (test<-0.49999995d*unit) { // singularity at south pole
					v.y = -2d * Mathd.Atan2 (y, x) * Mathd.Rad2Deg;
					v.x = -90d;
					v.z = 0d;
					return v;
				}
				v.y = Mathd.Rad2Deg * Mathd.Atan2 (2 * w * y + 2 * z * x, 1 - 2 * (x * x + y * y));     // Yaw
				v.x = Mathd.Rad2Deg * Mathd.Asin (2 * (w * x - y * z));                             // Pitch
				v.z = Mathd.Rad2Deg * Mathd.Atan2 (2 * w * z + 2 * x * y, 1 - 2 * (z * z + x * x));      // Roll

				for(int i=0; i<3; ++i){
					if(v[i] < 0d){
						v[i] += 360d;
					}
				}
//				v.y = Mathd.Rad2Deg * Mathd.Atan2 (2 * w * y + 2 * z * x, -sqw - sqz + sqx + sqy);     // Yaw
//				v.x = Mathd.Rad2Deg * Mathd.Asin (2 * (w * x - y * z));                             // Pitch
//				v.z = Mathd.Rad2Deg * Mathd.Atan2 (2 * w * z + 2 * x * y, sqw - sqz - sqx + sqy);      // Roll

//				Quaterniond q = new Quaterniond (w, z, x, y);
//				v.y = Mathd.Rad2Deg * Mathd.Atan2 (2 * q.x * q.w + 2 * q.y * q.z, 1 - 2 * (q.z * q.z + q.w * q.w));     // Yaw
//				v.x = Mathd.Rad2Deg * Mathd.Asin (2 * (q.x * q.z - q.w * q.y));                             // Pitch
//				v.z = Mathd.Rad2Deg * Mathd.Atan2 (2 * q.x * q.y + 2 * q.z * q.w, 1 - 2 * (q.y * q.y + q.z * q.z));      // Roll

				return v;

				
//				double sqx = x * x;
//				double sqy = y * y;
//				double sqz = z * z;
//				double sqw = w * w;
//				
//				Vector3d result = new Vector3d();
//				result.z = Mathd.Atan2(2.0d * (x*y + z*w),(sqx - sqy - sqz + sqw)) * Mathd.Rad2Deg;
//				result.x = Mathd.Atan2(2.0d * (y*z + x*w),(-sqx - sqy + sqz + sqw)) * Mathd.Rad2Deg;
//				result.y = Mathd.Asin(-2.0d * (x*z - y*w)) * Mathd.Rad2Deg;
//				
//				return result;
			}
		}
		
		public Vector3d eulerAnglesWrong {
			get {
				//				heading = atan2(2*qy*qw-2*qx*qz , 1 - 2*qy2 - 2*qz2)
				//					attitude = asin(2*qx*qy + 2*qz*qw)
				//						bank = atan2(2*qx*qw-2*qy*qz , 1 - 2*qx2 - 2*qz2)
				//						
				//						except when qx*qy + qz*qw = 0.5 (north pole)
				//						which gives:
				//						heading = 2 * atan2(x,w)
				//						bank = 0
				//						and when qx*qy + qz*qw = -0.5 (south pole)
				//						which gives:
				//						heading = -2 * atan2(x,w)
				//						bank = 0
				double heading, attitude, bank;
				attitude = Mathd.Asin (2*x*y + 2*z*w);
				double poleCheck = x * y + z * w;
				if(poleCheck == 0.5d){ 
					//	"North Pole"
					bank = 0;
					heading = 2 * Mathd.Atan2 (x, w);
				} else if(poleCheck == -0.5d){ 
					//	"South Pole"
					bank = 0;
					heading = -2 * Mathd.Atan2 (x, w);
				} else {
					heading = Mathd.Atan2 (2*y*w - 2*x*z, 1 - 2*y*y - 2*z*z);
					bank = Mathd.Atan2(2*x*w - 2*y*z,1 - 2*x*x - 2*z*z);
				}
				return new Vector3d(heading * Mathd.Rad2Deg, 
				                    attitude * Mathd.Rad2Deg, 
				                    bank * Mathd.Rad2Deg);
			}
		}

		public static Quaterniond FromToRotation(Vector3d from, Vector3d to){
			Quaterniond result = new Quaterniond();
			result.SetFromToRotation(from, to);
			return result;
		}
		public static Quaterniond FromToRotation(Quaterniond from, Quaterniond to){
			return Quaterniond.Inverse(from) * to;
		}
		public static Quaterniond Inverse(Quaternion rotation){
			Quaterniond result = new Quaterniond();
			result.x = -rotation.x;
			result.y = -rotation.y;
			result.z = -rotation.z;
			result.w = rotation.w;
			return result;
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
			from = from.normalized;
			to = to.normalized;
			double multA = 1 - t;
			double multB = t;

			if(Quaterniond.Dot(from, to) < 0d){
				from = -from;
			}
			
			return new Quaterniond(
				from.x * multA + to.x * multB,
				from.y * multA + to.y * multB,
				from.z * multA + to.z * multB,
				from.w * multA + to.w * multB
				).normalized;
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
			double angle = Quaterniond.Angle(from, to);
			double t = angle > maxDegreesDelta  ? maxDegreesDelta / angle : 1d;
//			Quaterniond difference = to.inverse * from;
			return Quaterniond.Slerp(from, to, t);
		}
		public static Quaterniond Slerp(Quaterniond from, Quaterniond to, double t){
			from = from.normalized;
			to = to.normalized;

			double dot = Quaterniond.Dot(from, to);
			if(dot < 0d){
				from = -from;
				dot = -dot;
			}
			if(dot > 0.9995){
				return Lerp (from, to, t);
			}
//			dot = Mathd.Clamp(dot, -1d, 1d);

			double theta =  Mathd.Acos(dot) * t;
			;
			double multA = Mathd.Cos(theta);
			double multB = Mathd.Sin(theta);

			Quaterniond c = to - from * dot;
			c.Normalize();

			return new Quaterniond(
				from.x * multA + c.x * multB,
				from.y * multA + c.y * multB,
				from.z * multA + c.z * multB,
				from.w * multA + c.w * multB
				).normalized;
		}

		public static Quaterniond operator -(Quaterniond q){
			return new Quaterniond(-q.x, -q.y, -q.z, -q.w);
		}

		public static bool operator !=(Quaterniond lhs, Quaterniond rhs){
			return !(lhs == rhs);
		}
		public static bool operator ==(Quaterniond lhs, Quaterniond rhs){
			return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w;
		}

		public static Quaterniond operator *(Quaterniond lhs, Quaterniond rhs){
			Quaterniond result = new Quaterniond();
			result.w = lhs.w * rhs.w - lhs.x*rhs.x - lhs.y*rhs.y - lhs.z*rhs.z;
			result.x = lhs.w * rhs.x + lhs.x*rhs.w + lhs.y*rhs.z - lhs.z*rhs.y;
			result.y = lhs.w * rhs.y - lhs.x*rhs.z + lhs.y*rhs.w + lhs.z*rhs.x;
			result.z = lhs.w * rhs.z + lhs.x*rhs.y - lhs.y*rhs.x + lhs.z*rhs.w;
			return result;
		}

		public static Quaterniond operator *(Quaterniond lhs, double rhs){
			Quaterniond result = new Quaterniond();
			result.w = lhs.w * rhs;
			result.x = lhs.x * rhs;
			result.y = lhs.y * rhs;
			result.z = lhs.z * rhs;
			return result;
		}
		
		public static Quaterniond operator *(Quaterniond lhs, Quaternion rhs){
			return lhs * new Quaterniond(rhs);
		}

		public static Quaterniond operator +(Quaterniond lhs, Quaterniond rhs){
			return new Quaterniond(lhs.x + rhs.x,
			                       lhs.y + rhs.y, 
			                       lhs.z + rhs.z, 
			                       lhs.w + rhs.w);
		}	
		public static Quaterniond operator -(Quaterniond lhs, Quaterniond rhs){
			return new Quaterniond(lhs.x - rhs.x,
			                       lhs.y - rhs.y, 
			                       lhs.z - rhs.z, 
			                       lhs.w - rhs.w);
		}

		public static Vector3d operator *(Quaterniond lhs, Vector3 rhs){
			return lhs.Multiply(rhs);
		}
		public static Vector3d operator *(Quaterniond lhs, Vector3d rhs){
			return lhs.Multiply(rhs);
		}
		public Vector3d Multiply(Vector3 rhs){
			return Multiply ((Vector3d)rhs);
		}
		public Vector3d Multiply(Vector3d rhs){
			double xx = x*x*2;
			double xy = x*y*2;
			double xz = x*z*2;
			double xw = x*w*2;
			double yy = y*y*2;
			double yz = y*z*2;
			double yw = y*w*2;
			double zz = z*z*2;
			double zw = z*w*2;
			return new Vector3d(
				//	ABij = A
				rhs.x * (1 - yy - zz) 	+ rhs.y * (xy - zw) 	+ rhs.z * (xz + yw),
				rhs.x * (xy + zw) 		+ rhs.y * (1 - xx - zz) + rhs.z * (yz - xw),
				rhs.x * (xz - yw) 		+ rhs.y * (yz + xw) 	+ rhs.z * (1 - xx - yy)
				);
		}
		public Vector3d InverseMultiply(Vector3 rhs){
			return InverseMultiply((Vector3d)rhs);
		}
		public Vector3d InverseMultiply(Vector3d rhs){
			double xx = x*x*2;
			double xy = x*y*2;
			double xz = x*z*2;
			double xw = -x*w*2;
			double yy = y*y*2;
			double yz = y*z*2;
			double yw = -y*w*2;
			double zz = z*z*2;
			double zw = -z*w*2;
			return new Vector3d(
				//	ABij = A
				rhs.x * (1 - yy - zz) 	+ rhs.y * (xy - zw) 	+ rhs.z * (xz + yw),
				rhs.x * (xy + zw) 		+ rhs.y * (1 - xx - zz) + rhs.z * (yz - xw),
				rhs.x * (xz - yw) 		+ rhs.y * (yz + xw) 	+ rhs.z * (1 - xx - yy)
				);
		}
//		Vector3d MultiplyOld(Vector3d rhs){
//			double[,] rotationMatrix = new double[,] {
//				{1 - 2*y*y - 2*z*z,	2*x*y - 2*z*w,	2*x*z + 2*y*w},
//				{2*x*y + 2*z*w, 1 - 2*x*x - 2*z*z,	2*y*z - 2*x*w},
//				{2*x*z - 2*y*w, 2*y*z + 2*x*w,	1 - 2*x*x - 2*y*y}};
//			Vector3d result = new Vector3d(
//				//	ABij = A
//				rhs.x * rotationMatrix[0,0] + rhs.y * rotationMatrix[0,1] + rhs.z * rotationMatrix[0,2],
//				rhs.x * rotationMatrix[1,0] + rhs.y * rotationMatrix[1,1] + rhs.z * rotationMatrix[1,2],
//				rhs.x * rotationMatrix[2,0] + rhs.y * rotationMatrix[2,1] + rhs.z * rotationMatrix[2,2]
//				);
//			return result;
//		}


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
			double sum = this.w * this.w + this.x * this.x + this.y * this.y + this.z * this.z;
//			double magnitudeInverse = Mathd.Sign (this.w) / System.Math.Sqrt(sum);
			double magnitudeInverse = 1d / System.Math.Sqrt(sum);
			this.w *= magnitudeInverse;
			this.x *= magnitudeInverse;
			this.y *= magnitudeInverse;
			this.z *= magnitudeInverse;
		}
//		public void NormalizeOld(){
//			double sum = 0;
//			for (int i = 0; i < 4; ++i){
//				sum += this[i] * this[i];
//			}
//			double magnitudeInverse = Mathd.Sign (this.w) / System.Math.Sqrt(sum);
//			for (int i = 0; i < 4; ++i){
//				this[i] *= magnitudeInverse;
//			}
//		}
		public void NormalizeSign(){
			if(this.w < 0d){
				this.w = -this.w;
				this.x = -this.x;
				this.y = -this.y;
				this.z = -this.z;
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

		public Quaternion ToFloatQuaternion(){
			return new Quaternion((float)x,(float)y,(float)z,(float)w);
		}

		/**
		 *	Tests if the values of the vector are approximately equal. The value of precisionBase is a multiple of the
		 *	base precision of a double, which is itself multiplied by the values to determine the tolerance. Any value
		 *	less than zero means "must be exactly zero".
		 */
		public static bool Approximately(Quaterniond lhs, Quaterniond rhs, double toleranceMultiplier){
			return Mathd.ApproximatelyRelative(lhs.x, rhs.x, toleranceMultiplier) 
				&& Mathd.ApproximatelyRelative(lhs.y, rhs.y, toleranceMultiplier)
					&& Mathd.ApproximatelyRelative(lhs.z, rhs.z, toleranceMultiplier)
			&& Mathd.ApproximatelyRelative(lhs.w, rhs.w, toleranceMultiplier);
		}
		
		public static bool ApproximatelyFixed(Quaterniond lhs, Quaterniond rhs, double toleranceBase){
			return Mathd.ApproximatelyFixed(lhs.x, rhs.x, toleranceBase) 
					&& Mathd.ApproximatelyFixed(lhs.y, rhs.y, toleranceBase)
					&& Mathd.ApproximatelyFixed(lhs.z, rhs.z, toleranceBase)
					&& Mathd.ApproximatelyFixed(lhs.w, rhs.w, toleranceBase);
		}

		static public implicit operator Quaterniond(Quaternion value) {
			return new Quaterniond(value);
		}
		static public explicit operator Quaternion(Quaterniond value)
		{
			return value.ToFloatQuaternion();
		}
	}
}

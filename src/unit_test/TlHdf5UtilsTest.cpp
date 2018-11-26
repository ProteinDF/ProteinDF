#include <string>
#include <vector>
#include "gtest/gtest.h"

#include "TlHdf5Utils.h"

static const std::string h5_path = "temp.h5";

TEST(TlHdf5Utils, createGroup) {
  TlHdf5Utils h5(h5_path);

  h5.createGroup("grp1");
  h5.createGroup("/grp2");

  EXPECT_TRUE(h5.hasGroup("grp1"));
  EXPECT_TRUE(h5.hasGroup("/grp1"));
  EXPECT_TRUE(h5.hasGroup("/grp1/"));
  EXPECT_TRUE(h5.hasGroup("grp1/"));
  EXPECT_TRUE(h5.hasGroup("grp2"));
  EXPECT_TRUE(h5.hasGroup("/grp2"));
  EXPECT_FALSE(h5.hasGroup("no_grp"));
  EXPECT_FALSE(h5.hasGroup("/no_grp"));
}

TEST(TlHdf5Utils, createGroup_withParents) {
  TlHdf5Utils h5(h5_path);

  h5.createGroup("/p1/grp1");
  h5.createGroup("/p2/p3/grp4");

  EXPECT_TRUE(h5.hasGroup("p1"));
  EXPECT_TRUE(h5.hasGroup("/p1"));
  EXPECT_TRUE(h5.hasGroup("/p1/grp1"));
  EXPECT_TRUE(h5.hasGroup("p1/grp1"));
  EXPECT_FALSE(h5.hasGroup("/p1/grp2"));
  EXPECT_FALSE(h5.hasGroup("/p2/grp1"));

  EXPECT_TRUE(h5.hasGroup("p2"));
  EXPECT_TRUE(h5.hasGroup("/p2"));
  EXPECT_TRUE(h5.hasGroup("/p2/p3"));
  EXPECT_TRUE(h5.hasGroup("p2/p3"));
  EXPECT_FALSE(h5.hasGroup("/p1/grp2"));
}

TEST(TlHdf5Utils, re_open) {
  TlHdf5Utils h5(h5_path);

  EXPECT_TRUE(h5.hasGroup("grp1"));
  EXPECT_FALSE(h5.hasGroup("grp_xxx"));
}

TEST(TlHdf5Utils, createDataSet_int) {
  TlHdf5Utils h5(h5_path);
  h5.createDataSet_int("large_empty_dataset_int", 1000);
}

TEST(TlHdf5Utils, createDataSet_double) {
  TlHdf5Utils h5(h5_path);
  h5.createDataSet_double("large_empty_dataset_double", 1000);
}

TEST(TlHdf5Utils, set_dataset_char) {
  TlHdf5Utils h5(h5_path);
  h5.write("group1/char1", 12);

  EXPECT_TRUE(h5.hasDataSet("/group1/char1"));
}

TEST(TlHdf5Utils, get_dataset_char) {
  TlHdf5Utils h5(h5_path);

  char char1;
  h5.get("/group1/char1", &char1);
  EXPECT_EQ(12, char1);
}

TEST(TlHdf5Utils, getset_int) {
  TlHdf5Utils h5(h5_path);
  h5.write("group1/int1", 123);

  EXPECT_TRUE(h5.hasDataSet("/group1/int1"));

  int int1;
  h5.get("/group1/int1", &int1);
  EXPECT_EQ(123, int1);
}

TEST(TlHdf5Utils, getset_long) {
  TlHdf5Utils h5(h5_path);
  h5.write("group1/long1", 12345678);

  EXPECT_TRUE(h5.hasDataSet("/group1/long1"));

  long long1;
  h5.get("/group1/long1", &long1);
  EXPECT_EQ(12345678, long1);
}

TEST(TlHdf5Utils, getset_int2long) {
  TlHdf5Utils h5(h5_path);
  h5.write("group1/int2", 555);

  EXPECT_TRUE(h5.hasDataSet("/group1/int2"));

  long a;
  h5.get("/group1/int2", &a);
  EXPECT_EQ(555, a);
}

TEST(TlHdf5Utils, getset_float) {
  TlHdf5Utils h5(h5_path);
  h5.write("grp1/float", 1.23);

  EXPECT_TRUE(h5.hasDataSet("/grp1/float"));

  float v;
  h5.get("/grp1/float", &v);
  EXPECT_FLOAT_EQ(1.23, v);
}

TEST(TlHdf5Utils, getset_double) {
  TlHdf5Utils h5(h5_path);
  h5.write("grp1/double", 3.14);

  EXPECT_TRUE(h5.hasDataSet("/grp1/double"));

  double v;
  h5.get("/grp1/double", &v);
  EXPECT_DOUBLE_EQ(3.14, v);
}

TEST(TlHdf5Utils, getset_int2double) {
  TlHdf5Utils h5(h5_path);
  h5.write("grp1/int2double", 777);

  EXPECT_TRUE(h5.hasDataSet("/grp1/int2double"));

  double v;
  h5.get("/grp1/int2double", &v);
  EXPECT_DOUBLE_EQ(777.0, v);
}

TEST(TlHdf5Utils, getset_string) {
  TlHdf5Utils h5(h5_path);
  h5.write("grp1/str", "hoge");

  EXPECT_TRUE(h5.hasDataSet("/grp1/str"));

  std::string string1;
  h5.get("/grp1/str", &string1);
  EXPECT_EQ("hoge", string1);
}

TEST(TlHdf5Utils, getset_strings) {
  TlHdf5Utils h5(h5_path);

  std::vector<std::string> strs(5);
  strs[0] = "This is a pen.";
  strs[1] = "Hello World!";
  strs[2] = "012345678901234567890123456789";
  strs[3] = "";
  strs[4] = "ProteinDF";

  h5.write("grp1/strs", strs);
  EXPECT_TRUE(h5.hasDataSet("/grp1/strs"));

  std::vector<std::string> out_strs;
  h5.get("/grp1/strs", &out_strs);
  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(strs[i], out_strs[i]);
  }
}

TEST(TlHdf5Utils, set_dataset_with_parents) {
  TlHdf5Utils h5(h5_path);
  h5.write("/p1/p2/ds1", "dataset1");

  EXPECT_TRUE(h5.hasGroup("/p1"));
  EXPECT_TRUE(h5.hasGroup("/p1/p2"));
  EXPECT_TRUE(h5.hasDataSet("/p1/p2/ds1"));
}

TEST(TlHdf5Utils, set_buffer_double) {
  const int size = 10;
  std::vector<double> vec_in(size);
  for (int i = 0; i < size; ++i) {
    vec_in[i] = 0.1 * i;
  }

  TlHdf5Utils h5(h5_path);
  h5.write("/buf/double", vec_in);
}

TEST(TlHdf5Utils, get_buffer_double) {
  int size = 5;  // <= 10
  TlHdf5Utils h5(h5_path);

  EXPECT_TRUE(h5.hasDataSet("/buf/double"));

  double* vec_out = new double[size];
  h5.get("/buf/double", vec_out, size);

  for (int i = 0; i < size; ++i) {
    EXPECT_EQ(0.1 * i, vec_out[i]);
  }

  delete[] vec_out;
  vec_out = NULL;
}

TEST(TlHdf5Utils, getset_vector_int) {
  const int size = 10;
  std::vector<int> vec_in(size);
  for (int i = 0; i < size; ++i) {
    vec_in[i] = i;
  }

  TlHdf5Utils h5(h5_path);
  h5.write("/vct/int", vec_in);

  EXPECT_TRUE(h5.hasDataSet("/vct/int"));

  std::vector<int> vec_out;
  h5.get("/vct/int", &vec_out);

  EXPECT_EQ(size, static_cast<int>(vec_out.size()));
  for (int i = 0; i < size; ++i) {
    EXPECT_EQ(vec_in[i], vec_out[i]);
  }
}

TEST(TlHdf5Utils, getset_vector_long) {
  int size = 10;
  std::vector<long> vec_in(size);
  for (int i = 0; i < size; ++i) {
    vec_in[i] = i;
  }

  TlHdf5Utils h5(h5_path);
  h5.write("/vct/long", vec_in);

  EXPECT_TRUE(h5.hasDataSet("/vct/long"));

  std::vector<long> vec_out;
  h5.get("/vct/long", &vec_out);

  EXPECT_EQ(size, static_cast<int>(vec_out.size()));
  for (int i = 0; i < size; ++i) {
    EXPECT_EQ(vec_in[i], vec_out[i]);
  }
}

TEST(TlHdf5Utils, getset_vector_float) {
  int size = 10;
  std::vector<float> vec_in(size);
  for (int i = 0; i < size; ++i) {
    vec_in[i] = i;
  }

  TlHdf5Utils h5(h5_path);
  h5.write("/vct/float", vec_in);

  EXPECT_TRUE(h5.hasDataSet("/vct/float"));

  std::vector<float> vec_out;
  h5.get("/vct/float", &vec_out);

  EXPECT_EQ(size, static_cast<int>(vec_out.size()));
  for (int i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(vec_in[i], vec_out[i]);
  }
}

TEST(TlHdf5Utils, getset_vector_double) {
  int size = 10;
  std::vector<double> vec_in(size);
  for (int i = 0; i < size; ++i) {
    vec_in[i] = i;
  }

  TlHdf5Utils h5(h5_path);
  h5.write("/vct/double", vec_in);

  EXPECT_TRUE(h5.hasDataSet("/vct/double"));

  std::vector<double> vec_out;

  h5.get("/vct/double", &vec_out);
  EXPECT_EQ(size, static_cast<int>(vec_out.size()));
  for (int i = 0; i < size; ++i) {
    EXPECT_DOUBLE_EQ(vec_in[i], vec_out[i]);
  }
}

TEST(TlHdf5Utils, getset_attr_int) {
  TlHdf5Utils h5(h5_path);
  h5.write("attr/int", "hogehoge");
  EXPECT_TRUE(h5.hasDataSet("/attr/int"));

  h5.setAttr("/attr/int", "last_modified", 2017);

  int buf;
  h5.getAttr("/attr/int", "last_modified", &buf);
  EXPECT_EQ(2017, buf);
}

TEST(TlHdf5Utils, getset_attr_long) {
  TlHdf5Utils h5(h5_path);
  h5.write("/attr/long", "hogehoge");
  EXPECT_TRUE(h5.hasDataSet("/attr/long"));

  h5.setAttr("/attr/long", "length", 1234567890);

  long buf;
  h5.getAttr("/attr/long", "length", &buf);
  EXPECT_EQ(1234567890, buf);
}

TEST(TlHdf5Utils, getset_attr_float) {
  TlHdf5Utils h5(h5_path);
  h5.write("/attr/float", "hogehoge");
  EXPECT_TRUE(h5.hasDataSet("/attr/float"));

  h5.setAttr("/attr/float", "param", 3.14);

  float buf;
  h5.getAttr("/attr/float", "param", &buf);
  EXPECT_FLOAT_EQ(3.14, buf);
}

TEST(TlHdf5Utils, getset_attr_double) {
  TlHdf5Utils h5(h5_path);
  h5.write("/attr/double", "hogehoge");
  EXPECT_TRUE(h5.hasDataSet("/attr/double"));

  h5.setAttr("/attr/double", "eff", 0.12345);

  double buf;
  h5.getAttr("/attr/double", "eff", &buf);
  EXPECT_DOUBLE_EQ(0.12345, buf);
}

TEST(TlHdf5Utils, get_selectedVector_int) {
  TlHdf5Utils h5(h5_path);

  const int size = 100;
  std::vector<int> v(size);
  for (int i = 0; i < size; ++i) {
    v[i] = i * 2;
  }
  h5.write("selvec/int", v);

  std::vector<unsigned long long> sel(5);
  sel[0] = 10;
  sel[1] = 17;
  sel[2] = 51;
  sel[3] = 3;
  sel[4] = 76;

  std::vector<int> out(sel.size());
  h5.getSelectedElements("/selvec/int", sel, &out);

  EXPECT_EQ(20, out[0]);
  EXPECT_EQ(34, out[1]);
  EXPECT_EQ(102, out[2]);
  EXPECT_EQ(6, out[3]);
  EXPECT_EQ(152, out[4]);
}

TEST(TlHdf5Utils, set_selectedVector_int) {
  TlHdf5Utils h5(h5_path);

  const int size = 100;
  std::vector<int> v(size);
  for (int i = 0; i < size; ++i) {
    v[i] = i * 2;
  }
  h5.write("selvec/int", v);

  std::vector<unsigned long long> sel(5);
  std::vector<int> data(sel.size());
  sel[0] = 10;
  sel[1] = 17;
  sel[2] = 51;
  sel[3] = 3;
  sel[4] = 76;
  data[0] = 0;
  data[1] = 1;
  data[2] = 2;
  data[3] = 3;
  data[4] = 4;
  h5.setSelectedElements("/selvec/int", sel, data);

  std::vector<int> out;
  h5.get("/selvec/int", &out);

  EXPECT_EQ(0, out[0]);
  EXPECT_EQ(2, out[1]);
  EXPECT_EQ(4, out[2]);
  EXPECT_EQ(0, out[10]);
  EXPECT_EQ(1, out[17]);
  EXPECT_EQ(2, out[51]);
  EXPECT_EQ(3, out[3]);
  EXPECT_EQ(4, out[76]);
}

TEST(TlHdf5Utils, get_selected_double) {
  TlHdf5Utils h5(h5_path);
  // h5.createGroup("/get_selected");

  const int size = 100;
  std::vector<double> v(size);
  for (int i = 0; i < size; ++i) {
    v[i] = double(i) * 2.0;
  }
  h5.write("/get_selected/double", v);

  std::vector<unsigned long long> sel(5);
  sel[0] = 10;
  sel[1] = 17;
  sel[2] = 51;
  sel[3] = 3;
  sel[4] = 76;

  std::vector<double> out(sel.size());
  h5.getSelectedElements("/get_selected/double", sel, &out);

  EXPECT_DOUBLE_EQ(20.0, out[0]);
  EXPECT_DOUBLE_EQ(34.0, out[1]);
  EXPECT_DOUBLE_EQ(102.0, out[2]);
  EXPECT_DOUBLE_EQ(6.0, out[3]);
  EXPECT_DOUBLE_EQ(152.0, out[4]);
}

TEST(TlHdf5Utils, set_selected_double) {
  TlHdf5Utils h5(h5_path);
  h5.createGroup("/set_selected");

  const int size = 100;
  std::vector<double> v(size);
  for (int i = 0; i < size; ++i) {
    v[i] = double(i) * 2.0;
  }
  h5.write("/set_selected/double", v);

  std::vector<unsigned long long> sel(5);
  std::vector<double> data(sel.size());
  sel[0] = 10;
  sel[1] = 17;
  sel[2] = 51;
  sel[3] = 3;
  sel[4] = 76;
  data[0] = 0.0;
  data[1] = 1.0;
  data[2] = 2.0;
  data[3] = 3.0;
  data[4] = 4.0;
  h5.setSelectedElements("/set_selected/double", sel, data);

  std::vector<double> out;
  h5.get("/set_selected/double", &out);

  EXPECT_DOUBLE_EQ(0.0, out[0]);
  EXPECT_DOUBLE_EQ(2.0, out[1]);
  EXPECT_DOUBLE_EQ(4.0, out[2]);
  EXPECT_DOUBLE_EQ(0.0, out[10]);
  EXPECT_DOUBLE_EQ(1.0, out[17]);
  EXPECT_DOUBLE_EQ(2.0, out[51]);
  EXPECT_DOUBLE_EQ(3.0, out[3]);
  EXPECT_DOUBLE_EQ(4.0, out[76]);
}

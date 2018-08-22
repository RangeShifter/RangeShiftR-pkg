object frmMain: TfrmMain
  Left = 0
  Top = 0
  Width = 1308
  Height = 784
  HorzScrollBar.Range = 1308
  VertScrollBar.Range = 1200
  Anchors = [akLeft, akBottom]
  Caption = 'RangeShifter_v2.0'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  Menu = MainMenu1
  OldCreateOrder = False
  Visible = True
  PixelsPerInch = 96
  TextHeight = 13
  object PopLegend: TImage
    Left = 634
    Top = 63
    Width = 10
    Height = 300
  end
  object LabelPop: TLabel
    Left = 624
    Top = 25
    Width = 72
    Height = 32
    Caption = 'Population Size'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    Visible = False
    WordWrap = True
  end
  object LabelMinPop: TLabel
    Left = 650
    Top = 64
    Width = 7
    Height = 16
    BiDiMode = bdLeftToRight
    Caption = '0'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentBiDiMode = False
    ParentFont = False
    ParentShowHint = False
    ShowHint = True
    Visible = False
    WordWrap = True
  end
  object LabelMaxPop: TLabel
    Left = 650
    Top = 348
    Width = 28
    Height = 16
    BiDiMode = bdLeftToRight
    Caption = 'max '
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentBiDiMode = False
    ParentFont = False
    ParentShowHint = False
    ShowHint = True
    Visible = False
  end
  object LandScrollBox: TScrollBox
    Left = 0
    Top = 0
    Width = 622
    Height = 650
    TabOrder = 4
    object LandImage: TImage
      Left = 0
      Top = 0
      Width = 600
      Height = 1000
    end
    object CommImage: TImage
      Left = 0
      Top = 0
      Width = 600
      Height = 1000
    end
    object SpDistImage: TImage
      Left = 0
      Top = 0
      Width = 600
      Height = 1000
      Transparent = True
    end
    object MovtPaintBox: TPaintBox
      Left = 0
      Top = 0
      Width = 600
      Height = 1000
    end
  end
  object ChartOccSuit: TChart
    Left = 700
    Top = 575
    Width = 561
    Height = 180
    BackWall.Visible = False
    Legend.Alignment = laTop
    Legend.Font.Height = -13
    Legend.LegendStyle = lsSeries
    Legend.Shadow.Visible = False
    Legend.TopPos = 47
    Legend.Transparent = True
    Legend.Visible = False
    MarginBottom = 0
    MarginLeft = 0
    MarginRight = 5
    MarginTop = 5
    RightWall.Visible = True
    Title.Font.Color = clBlack
    Title.Font.Height = -13
    Title.Font.Style = [fsBold]
    Title.Text.Strings = (
      'Proportion of suitable cells occupied')
    BottomAxis.Title.Caption = 'Years'
    BottomAxis.Title.Font.Height = -13
    BottomAxis.Title.Font.Style = [fsBold]
    DepthAxis.Automatic = False
    DepthAxis.AutomaticMaximum = False
    DepthAxis.AutomaticMinimum = False
    DepthAxis.Maximum = -0.100000000000000600
    DepthAxis.Minimum = -1.100000000000001000
    DepthTopAxis.Automatic = False
    DepthTopAxis.AutomaticMaximum = False
    DepthTopAxis.AutomaticMinimum = False
    DepthTopAxis.Maximum = -0.100000000000000600
    DepthTopAxis.Minimum = -1.100000000000001000
    LeftAxis.ExactDateTime = False
    LeftAxis.Grid.SmallDots = True
    LeftAxis.Title.Font.Height = -13
    LeftAxis.Title.Font.Style = [fsBold]
    RightAxis.Automatic = False
    RightAxis.AutomaticMaximum = False
    RightAxis.AutomaticMinimum = False
    RightAxis.Grid.Visible = False
    RightAxis.Title.Font.Height = -13
    RightAxis.Title.Font.Style = [fsBold]
    RightAxis.Visible = False
    TopAxis.Automatic = False
    TopAxis.AutomaticMaximum = False
    TopAxis.AutomaticMinimum = False
    TopAxis.Visible = False
    View3D = False
    View3DWalls = False
    BevelInner = bvLowered
    TabOrder = 2
    Visible = False
    DefaultCanvas = 'TGDIPlusCanvas'
    ColorPaletteIndex = 13
    object OccSuitmean: TLineSeries
      SeriesColor = clBlack
      Title = 'Occupied Cells / Suitable Cells'
      Brush.BackColor = clDefault
      LinePen.Width = 2
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Transparency = 10
    end
    object OccSuitplusSE: TLineSeries
      SeriesColor = clBlack
      Brush.BackColor = clDefault
      LinePen.Style = psDot
      LinePen.Width = 2
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object OccSuitminusSE: TLineSeries
      SeriesColor = clBlack
      Brush.BackColor = clDefault
      LinePen.Style = psDot
      LinePen.Width = 2
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
  end
  object Memo1: TMemo
    Left = 700
    Top = 8
    Width = 561
    Height = 121
    ReadOnly = True
    ScrollBars = ssVertical
    TabOrder = 1
  end
  object ChartPop: TChart
    Left = 700
    Top = 151
    Width = 561
    Height = 194
    Legend.Alignment = laTop
    Legend.Font.Height = -13
    Legend.LegendStyle = lsSeries
    Legend.TopPos = 47
    Legend.Transparent = True
    Legend.Visible = False
    MarginBottom = 0
    MarginLeft = 0
    MarginRight = 0
    MarginTop = 5
    Title.Text.Strings = (
      '')
    Title.Visible = False
    Title.AdjustFrame = False
    BottomAxis.Grid.Visible = False
    BottomAxis.Title.Caption = 'Years'
    BottomAxis.Title.Font.Height = -13
    BottomAxis.Title.Font.Style = [fsBold]
    DepthAxis.Automatic = False
    DepthAxis.AutomaticMaximum = False
    DepthAxis.AutomaticMinimum = False
    DepthAxis.Maximum = 0.500000000000000000
    DepthAxis.Minimum = -0.500000000000000000
    DepthTopAxis.Automatic = False
    DepthTopAxis.AutomaticMaximum = False
    DepthTopAxis.AutomaticMinimum = False
    DepthTopAxis.Maximum = 0.500000000000000000
    DepthTopAxis.Minimum = -0.500000000000000000
    LeftAxis.Title.Caption = 'Population size'
    LeftAxis.Title.Font.Color = clRed
    LeftAxis.Title.Font.Height = -13
    LeftAxis.Title.Font.Style = [fsBold]
    RightAxis.Grid.Visible = False
    RightAxis.Title.Caption = 'Occupied cells'
    RightAxis.Title.Font.Color = clBlue
    RightAxis.Title.Font.Height = -13
    RightAxis.Title.Font.Style = [fsBold]
    TopAxis.Automatic = False
    TopAxis.AutomaticMaximum = False
    TopAxis.AutomaticMinimum = False
    TopAxis.Visible = False
    View3D = False
    View3DWalls = False
    BevelInner = bvLowered
    TabStop = False
    TabOrder = 0
    DefaultCanvas = 'TGDIPlusCanvas'
    ColorPaletteIndex = 13
    object Population: TLineSeries
      SeriesColor = clRed
      Brush.BackColor = clDefault
      LinePen.Width = 2
      Pointer.Brush.Gradient.EndColor = 10708548
      Pointer.Gradient.EndColor = 10708548
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object Cells: TLineSeries
      SeriesColor = clBlue
      VertAxis = aRightAxis
      Brush.BackColor = clDefault
      LinePen.Width = 2
      Pointer.Brush.Gradient.EndColor = 3513587
      Pointer.Gradient.EndColor = 3513587
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
  end
  object ChartNoise: TChart
    Left = 700
    Top = 367
    Width = 561
    Height = 180
    Legend.Alignment = laTop
    Legend.Font.Height = -13
    Legend.LegendStyle = lsSeries
    Legend.TopPos = 47
    Legend.Transparent = True
    Legend.Visible = False
    MarginBottom = 0
    MarginLeft = 0
    MarginRight = 5
    MarginTop = 5
    RightWall.Visible = True
    Title.AutoSize = False
    Title.Font.Color = clBlack
    Title.Font.Height = -13
    Title.Font.Style = [fsBold]
    Title.Text.Strings = (
      'Environmental noise time series')
    Title.AdjustFrame = False
    Title.VertMargin = 0
    BottomAxis.Grid.Visible = False
    BottomAxis.Title.Caption = 'Years'
    BottomAxis.Title.Font.Height = -13
    BottomAxis.Title.Font.Style = [fsBold]
    DepthAxis.Automatic = False
    DepthAxis.AutomaticMaximum = False
    DepthAxis.AutomaticMinimum = False
    DepthAxis.Maximum = 0.500000000000000000
    DepthAxis.Minimum = -0.500000000000000000
    DepthTopAxis.Automatic = False
    DepthTopAxis.AutomaticMaximum = False
    DepthTopAxis.AutomaticMinimum = False
    DepthTopAxis.Maximum = 0.500000000000000000
    DepthTopAxis.Minimum = -0.500000000000000000
    LeftAxis.Title.Caption = 'epsilon'
    LeftAxis.Title.Font.Height = -13
    LeftAxis.Title.Font.Style = [fsBold]
    RightAxis.Grid.Visible = False
    RightAxis.Title.Caption = 'Occupied Cells'
    RightAxis.Title.Font.Height = -13
    RightAxis.Title.Font.Style = [fsBold]
    TopAxis.Automatic = False
    TopAxis.AutomaticMaximum = False
    TopAxis.AutomaticMinimum = False
    TopAxis.Visible = False
    View3D = False
    View3DWalls = False
    BevelInner = bvLowered
    TabStop = False
    TabOrder = 3
    Visible = False
    DefaultCanvas = 'TGDIPlusCanvas'
    ColorPaletteIndex = 13
    object EnvNoise: TLineSeries
      SeriesColor = clRed
      Brush.BackColor = clDefault
      LinePen.Width = 2
      Pointer.Brush.Gradient.EndColor = 10708548
      Pointer.Gradient.EndColor = 10708548
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
  end
  object BtnZoomIn: TButton
    Left = 16
    Top = 668
    Width = 62
    Height = 33
    Caption = 'Zoom in'
    Default = True
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ModalResult = 1
    ParentFont = False
    TabOrder = 5
    TabStop = False
    OnClick = BtnZoomInClick
  end
  object BtnZoomOut: TButton
    Left = 98
    Top = 668
    Width = 60
    Height = 33
    Caption = 'Zoom out'
    Default = True
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ModalResult = 1
    ParentFont = False
    TabOrder = 6
    TabStop = False
    OnClick = BtnZoomOutClick
  end
  object MainMenu1: TMainMenu
    Left = 88
    Top = 8
    object File1: TMenuItem
      Caption = 'File      '
      object SetDirectory: TMenuItem
        Caption = 'Set Directory'
        OnClick = SetDirectoryClick
      end
      object BatchMode: TMenuItem
        Caption = 'Batch Mode'
        Enabled = False
        OnClick = BatchModeClick
      end
    end
    object LandMenu: TMenuItem
      Caption = 'Landscape'
      Enabled = False
      OnClick = LandMenuClick
      object RasterLand: TMenuItem
        Caption = 'Import Raster'
        OnClick = RasterLandClick
      end
      object Artificial: TMenuItem
        Caption = 'Generate Artificial Landscape'
        OnClick = ArtificialClick
      end
      object EnvGradient: TMenuItem
        Caption = 'Environmental Gradient'
        Enabled = False
        OnClick = EnvGradientClick
      end
    end
    object Parameterset: TMenuItem
      Caption = 'Parameters setting'
      Enabled = False
      object SpeciesMenu: TMenuItem
        Caption = 'Species'
        Enabled = False
        OnClick = SpeciesMenuClick
      end
      object GeneticsMenu: TMenuItem
        Caption = 'Genetics'
        Enabled = False
        OnClick = GeneticsMenuClick
      end
      object Simulations: TMenuItem
        Caption = 'Simulations'
        Enabled = False
        OnClick = SimulationsClick
      end
    end
    object Run1: TMenuItem
      Bitmap.Data = {
        36030000424D3603000000000000360000002800000010000000100000000100
        18000000000000030000120B0000120B00000000000000000000FF00FFFF00FF
        FF00FFFF00FFFF00FF6D3327853C1395440D96450D873D12703425FF00FFFF00
        FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF70372A70372ACD7727E8AD70F3
        CCA1F4CDA3E9B176D07C2C6F35296F3529FF00FFFF00FFFF00FFFF00FFFF00FF
        86411DC0620BF0C292FFFEFAFDFAF6F5E3D1F5E2D0FDF8F4FFFFFDF2C99EC669
        117B3A21FF00FFFF00FFFF00FF8C451CC16107F7DBBDFFFEFEE0A46BCE6D13C7
        5C00C96100CE6E12DE9D5FFDFAF7FAE5CCC6680D6F3528FF00FFFF00FF8C451C
        ECBD8BFFFFFFDA8F43C65600FFFFFFDC9751C75B00CA6200C75B00D58333FDFA
        F8F3CB9F6F3528FF00FFA04D10CE7721FFFDFBE8B684D06B04D06B07FFFFFFFF
        FFFFE1A870C95E00CA6300C75B00DFA061FFFFFFCF7B28703525AF5507E5AA6F
        FFFFFFDD8F3FDA8128D87B1EFFFFFFFFFFFFFFFFFFEAC198CC6708C95F00CE6E
        0DFDFAF6E9B175703525BB5F0AF0CAA1FCF4EDE19343E29242DF8A34FFFFFFFF
        FFFFFFFFFFFFFFFFF3DEC6CF7017C95F00F5E3D0F3CEA4703525C1650FF2CDA6
        FDF7F0E9A158E9A056E6994AFFFFFFFFFFFFFFFFFFFFFFFFEBC39BCD6A0DC961
        00F6E6D4F3CCA1703525C1640DEEBC88FFFFFFF1B87CF0AE69EEA75FFFFFFFFF
        FFFFFFFFFFE7B17ACE6902C96100CF7111FEFCFAE7AC6D703525BF6006E5A059
        FFFDFAFBE0C4F8BA7BF4B471FFFFFFFFFFFFE8AB6DD87B1DD27414C85C00E2AA
        71FFFFFECC7520703525FF00FFC36204FAD9B8FFFFFFFEDCB8F7B877FFFFFFEB
        AD6DE08D37DA8228D06B05DB924AFFFFFFEFC08C6B342CFF00FFFF00FFC36204
        E79E55FEEBD7FFFFFFFBDFC3F1B578E89F55E29445DE9142E8B786FFFFFFF6D8
        B7BE5F066B342CFF00FFFF00FFFF00FFC6670CE69E55FAD9B6FFFBF6FFFFFFFE
        F8F2FDF6EFFFFFFFFEF9F2ECB884BE5F09753826FF00FFFF00FFFF00FFFF00FF
        FF00FFC06005C06005E49F5AEEBA86F2CAA0F0C599E4A768CC741E783A27783A
        27FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFB65C0AB86012B9
        6113B25A0FA24F0E8E451AFF00FFFF00FFFF00FFFF00FFFF00FF}
      Caption = 'Run'
      Enabled = False
      OnClick = Run1Click
    end
    object Pause: TMenuItem
      Bitmap.Data = {
        36030000424D3603000000000000360000002800000010000000100000000100
        18000000000000030000120B0000120B00000000000000000000FF00FFFF00FF
        FF00FFFF00FFFF00FF7D2C057D2C057D2C057D2C057D2C057D2C05FF00FFFF00
        FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF552117552117BF5C15E19855EF
        BD8AF0BF8CE29D5BC26218541F16541F16FF00FFFF00FFFF00FFFF00FFFF00FF
        6C290EAF4704EBB179FFFEF8FEF8F4F2DAC4F2D9C2FCF6F0FFFFFCEEBA87B64E
        07602311FF00FFFF00FFFF00FF732D0DB04602F4D0ABFFFFFED48846BC4B02B8
        4300B84300BB4901D17F3BFCFAF4F8DDBDB64D05541F16FF00FFFF00FF732D0D
        E6AB72FFFFFFCC7228CA6B23D89151C56012C45C0DD89152CB6E28C7661DFCF8
        F4EFBC88541F16FF00FFB04A06C05C11FFFCFAE1A36ABF4800FFFFFFFFFFFFD4
        8742CF7930FFFFFFFFFFFFB53900D58945FFFFFFC160167B2C04B04A06DD9554
        FFFFFFD37628CC600EFFFFFFFFFFFFD98E49D37F36FFFFFFFFFFFFB84100C053
        05FCF8F3E29D5A7B2C04B04A06EBBB8AFBF0E7D87A2BD77422FFFFFFFFFFFFDE
        9655D98942FFFFFFFFFFFFBA4400BB4600F2DAC2EFC08E7B2C04B04A06EEBF90
        FCF4EBE28A3EE18435FFFFFFFFFFFFE29E5FDD914DFFFFFFFFFFFFBA4400BB47
        00F3DEC7EFBD8A7B2C04B04A06E9AA6EFFFFFFEDA462EA9547FFFFFFFFFFFFE6
        A569E19856FFFFFFFFFFFFB84200C15506FEFBF8E097527B2C04B04A06DD893F
        FFFCF8FAD7B4F4A359FFFFFFFFFFFFEAAB72E59E5EFFFFFFFFFFFFB63C00D995
        56FFFFFEBD5A107B2C04FF00FFB24701F8CEA5FFFFFFFED0A3F7B77BF2B780E5
        934CDD8438DE9350CE6E21CF782FFFFFFFEAAF73501F18FF00FFFF00FFB24701
        E0873BFEE5CBFFFFFFFAD3ADEDA05AE1893CD9792AD17021E09F64FFFFFFF3CC
        A4AC4402501F18FF00FFFF00FFFF00FFB64C04DE873BF8CEA3FFFAF3FFFFFFFE
        F6EEFCF3EAFFFFFFFEF7EEE6A56AAC44035A2214FF00FFFF00FFFF00FFFF00FF
        FF00FFAF4501AF4501DC8840E9A76CEEBB89EBB581DC914DBD590F5D23155D23
        15FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFA64608A64608A6
        4608A64608A64608A64608FF00FFFF00FFFF00FFFF00FFFF00FF}
      Caption = 'Pause'
      Enabled = False
      OnClick = PauseClick
    end
    object Stop: TMenuItem
      Bitmap.Data = {
        36030000424D3603000000000000360000002800000010000000100000000100
        18000000000000030000120B0000120B00000000000000000000FF00FFFF00FF
        FF00FFFF00FFFF00FF000288010893010B99010C99010893000389FF00FFFF00
        FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF000186010D9D021CAF021FB402
        1FB5021FB5021FB2021CB0010F9F000287FF00FFFF00FFFF00FFFF00FFFF00FF
        00038A0118B2021FB7021EB1021DB1021DB1021DB1021DB1021EB2021FB40219
        AC00048EFF00FFFF00FFFF00FF0001830118BB0220CC011CBF0015B4011AB002
        1DB1021DB1011CB00015AD011BB0021FB4021AAC000287FF00FFFF00FF010CA7
        0121E0011CD30726CC4966D70B28BC0018B00019AF0622B44A66CE0D2BB7011B
        B0021FB5010F9FFF00FF000187011CDC0120ED0017DC3655E2FFFFFFA4B4ED05
        20BB0119B28B9EE1FFFFFF4E6ACF0014AC021EB2021CB000038900069A0120F8
        011FF6001DE9031FE1738BEEFFFFFFA0B1ED93A5E7FFFFFF91A4E20823B4011B
        B0021DB1021FB4010895020CAA0A2EFE0323FB011FF6001CEB0018E1788FF0FF
        FFFFFFFFFF96A7EA021BB50019AF021DB1021DB10220B5010C99040EAC294DFE
        0D30FB011FFA001CF7011CEE8EA1F4FFFFFFFFFFFFA6B6EE0520C20018B6021D
        B1021DB10220B5010B980208A04162FB2F51FC001EFA0725FA8AA0FEFFFFFF8E
        A3F67991F2FFFFFFA3B4EE0C29C6011BB8021DB4021FB2000793000189314FEF
        7690FF0F2DFA3354FBFFFFFF91A5FE021EF30017E7738BF3FFFFFF4765E00016
        C2021FBD021CB2000288FF00FF0C1BBC819AFF728BFE1134FA3456FB0526FA00
        1CFA001CF40220ED3353ED0625DA011DD00220CB010DA1FF00FFFF00FF000189
        2943E6A5B7FF849AFC2341FB0323FA011FFA011FFA001EF7011BEE021FE50121
        E20118BF000184FF00FFFF00FFFF00FF01038F2A45E693A9FFABBBFF758FFE49
        69FC3658FB3153FC2346FC092CF70118CB00038BFF00FFFF00FFFF00FFFF00FF
        FF00FF0001890F1DBF3E5BF36B87FE728CFF5E7BFE395BFB1231EB010FB50001
        84FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF000189030AA306
        11B2050FB10107A0000188FF00FFFF00FFFF00FFFF00FFFF00FF}
      Caption = 'Stop'
      Enabled = False
      OnClick = StopClick
    end
    object Refresh: TMenuItem
      Bitmap.Data = {
        36030000424D3603000000000000360000002800000010000000100000000100
        18000000000000030000120B0000120B00000000000000000000FF00FFFF00FF
        FF00FFFF00FFFF00FFA23F08FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00
        FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFA13E08B04D06A13E08FF
        00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF
        FF00FFA13E08BD5804BF5B03A13E08FF00FFFF00FFFF00FFFF00FFFF00FFFF00
        FFFF00FFFF00FFFF00FFFF00FFFF00FFA13D08B34F05C86401A13E08A13E08FF
        00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFA13D08
        AE4C09C96502A44108A13E08FF00FFA13E08A13E08A03D08A03D08A03D08A03D
        08A03D08A13E08A13E08FF00FFA44209D97E1CAA4809A13E08FF00FFA13E08B0
        4D06C66102C66002C56002C56003C45F03C45F03C56003A03E08A03C06CE7626
        C96F1FA03C07FF00FFFF00FFFF00FFA13E08A13E08A13E08A13E08AD4B06CC68
        00C96402D06B00A03D089F3B06ED9D40A13E08FF00FFFF00FFFF00FFFF00FFFF
        00FFFF00FFFF00FF9F3C08BD5904C66202AA4606D06B00A03D089E3A05FFBF55
        9D3905FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFA03D08B24E05CD6900A23F
        09A13E08D06B00A03D089E3B05FAAE4BB151149F3C06FF00FFFF00FFFF00FFFF
        00FFA03D08BC5703CE6900A94607A03D08A13E08D16C00A03D089F3C06D27C2D
        E6953C9D3905A13E08FF00FF9F3C06A34109CC6D13CB6806A64207A13E08FF00
        FFA13E08D16C00A03D08FF00FFA5430AEE9E41E7953CB85A18B95A19C56A21ED
        9A37CA6E1BA44109A13E08FF00FFFF00FFA13E08D16C00A03D08FF00FFA03D07
        9F3B06BF631EE08E38EE9E41DA8532B15213A03C07A13E08FF00FFFF00FFFF00
        FFA13E08C96401A13E08FF00FFFF00FFFF00FFA03C069F3C069F3B069F3C06A0
        3D07FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFA13E08FF00FFFF00FFFF00FF
        FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00
        FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF
        00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF}
      Caption = 'Refresh'
      Enabled = False
      OnClick = RefreshClick
    end
  end
  object OpenDialog1: TOpenDialog
    Left = 16
    Top = 8
  end
  object OpenDialog2: TOpenDialog
    Left = 16
    Top = 64
  end
end

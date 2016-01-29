define(
    ['jquery', 'app/ui/detail', 'app/ui/layout', 'app/util/utils', 'endless', 'bootstrap'],
    function($, Detail, Layout, Util) {

    var detail, $search, $togglePinned,
        $toggleCluster, $clusterSelection, selectedCluster,
        $toggleTimespan, $timespanSelection,
        loadThumbsXhr;

    function selectCluster(selection) {
        if (selection == null || selection == 'none') {
            selection = 'none';
        } else {
            selection = parseInt(selection);
        }

        $clusterSelection.removeClass('cluster-' + selectedCluster);
        selectedCluster = selection;
        if (selectedCluster == 'none') {
            $clusterSelection.hide();
            $toggleCluster.find('.cluster-icon').show();
        } else if (0 <= selection && selection <= 3) {
            $clusterSelection
                .css({display: 'inline-block'})
                .addClass('cluster-' + selectedCluster);
            $toggleCluster.find('.cluster-icon').hide();
        } else {
            throw 'Unexpected cluster selection';
        }
        $clusterSelection.data('value', selection);
    }

    function selectTimespan(selection) {
        var value = (function(s) {
            switch (s) {
                case 7 : return 'W';
                case 30 : return '1m';
                case 60 : return '2m';
            }
            throw 'Unexpected timespan selection';
        })(selection);

        $timespanSelection.html(value).data('value', selection);
    }

    function getControlsStates() {
        var cluster = $clusterSelection.data('value'),
            filters = [];

        cluster = cluster == 'none' ? null : parseInt(cluster);

        $('.filter-group').each(function() {
            var $group = $(this),
                values = [];
            $group.find('li a.active').each(function() {
                values.push($(this).data('id'));
            });
            if (values.length) {
                filters.push({id: $group.data('id'), pk: values});
            }
        });

        return {
            time_span: parseInt($timespanSelection.data('value')),
            cluster: cluster,
            pin: ($('#toggle-pinned').hasClass('active')),
            search_query: null, //$('#search-input').val(),
            filters: filters
        };
    }

    function updateControlsStates(data) {
        // Cluster
        if (data.hasOwnProperty('cluster')) {
            selectCluster(data['cluster'])
        }
        // Time span
        if (data.hasOwnProperty('time_span')) {
            selectTimespan(data['time_span']);
        }
        // Pinned
        if (data.hasOwnProperty('pin')) {
            Util.toggleClass($togglePinned, 'active', data['pin'])
        }
        // Search
        if (data.hasOwnProperty('search_query')) {
            $('#search-input').val(data['search_query']);
        }
        // Filters
        if (data.hasOwnProperty('filters')) {
            $(data['filters']).each(function (i, filter) {
                // TODO
            });
        }
    }

    function loadThumbs() {
        if (loadThumbsXhr) {
            loadThumbsXhr.abort();
        }

        loadThumbsXhr = $.ajax({
            url: '/feed/stream2/filter',
            data: getControlsStates(),
            dataType: 'xml',
            method: 'GET'
        });
        loadThumbsXhr.success(function(xml) {
            var data = $(xml).find('data');
            if (data.length) {
                updateControlsStates(JSON.parse(data.text()));
            }

            var list = $(xml).find('thumb-list');
            if (list.length) {
                $('#list .endless_data').html($(list.text()));
            }

            loadThumbsXhr = null;
        });
        loadThumbsXhr.error(function() {
            console.log('Fail to load thumbs.');
        });
    }

    $(function() {

        $search = $('#search');
        $togglePinned = $('#toggle-pinned');

        $toggleCluster = $('#toggle-cluster');
        $clusterSelection = $('#cluster-selection');

        $toggleTimespan = $('#toggle-timespan');
        $timespanSelection = $('#timespan-selection');

        // Toggle pinned
        $togglePinned.on('click', function(e) {
            Util.toggleClass($(e.delegateTarget), 'active');
            loadThumbs();
        });

        // Toggle search bar
        $('#toggle-search, #close-search').on('click', function() {
            Util.toggleClass($search, 'opened');
        });

        // Search bar active
        $('#search-input')
            .on('focus', function(e) {
                e.stopPropagation();
                $(e.delegateTarget).parents('form').addClass('active');
            })
            .on('blur', function(e) {
                e.stopPropagation();
                $(e.delegateTarget).parents('form').removeClass('active');
            });

        // Toggle selector
        $('.selector button').on('click', function(e) {
            //e.stopPropagation();
            var $selector = $(e.target).parents('.selector').eq(0);
            Util.toggleClass($selector, 'opened');
        });

        // Cluster selection
        $('#cluster').on('click', '.choices .cluster-value', function(e) {
            selectCluster($(e.target).closest('.cluster-value').data('cluster'));
            $toggleCluster.trigger('click');
            loadThumbs();
        });

        // Timespan selection
        $('#timespan').on('click', '.choices a', function(e) {
            selectTimespan($(e.target).closest('a').data('timespan'));
            $toggleTimespan.trigger('click');
            loadThumbs();
        });

        // Close selectors on click out
        $(window).on('click', function(e) {
            if (0 == $(e.target).closest('.selector').length) {
                $('.selector').removeClass('opened');
            }
        });


        // Filters
        $('#filter-flap')
            .on('click', '.filter-toggle', function(e) {
                var $group = $(e.target).parents('.filter-group').eq(0);
                if (Util.toggleClass($group, 'active')) {
                    $group.find('.filter-filters').collapse('show');
                } else {
                    $group.find('.filter-filters').collapse('hide');
                }
            })
            .on('click', '.filter-group ul a', function(e) {
                Util.toggleClass($(e.target), 'active');
                loadThumbs();
            });


        // Thumbs
        $('#list, #detail')
            .on('click', '.thumb-pin', function(e) {
                e.stopPropagation();

                var $button = $(this);
                var $thumb = $(e.target).parents('.thumb').eq(0);

                var pinXhr = $.ajax({
                    type: 'POST',
                    url: '/user/paper/pin',
                    data: {'pk': $thumb.data('id'), 'source': window.location.pathname}
                });
                pinXhr.success(function(json) {
                    if (json.hasOwnProperty('is_liked')) {
                        $button.toggleClass('active', json['is_liked']);
                    }
                });
                pinXhr.fail(function() {
                    console.log('ERROR');
                });

                return false;
            })
            .on('click', '.thumb-remove', function(e) {
                e.stopPropagation();

                var $thumb = $(e.target).parents('.thumb').eq(0);

                var pinXhr = $.ajax({
                    type: 'POST',
                    url: '/user/paper/ban',
                    data: {'pk': $thumb.data('id'), 'source': window.location.pathname}
                });
                pinXhr.success(function(json) {
                    if (json.hasOwnProperty('is_liked') && json['is_ticked']) {
                        $thumb.remove();
                    }
                });
                pinXhr.fail(function() {
                    console.log('ERROR');
                });

                return false;
            });

        // Endless scroll
        $.endlessPaginate({
            paginateOnScroll: true,
            paginateOnScrollMargin: 10
        });


        // Detail
        detail = new Detail();
        $(detail)
            .on('etalia.detail.loading', function() {
                Layout.setBusy();
            })
            .on('etalia.detail.loaded', function() {
                Layout.setAvailable();
            });
        $('.document').on('click', '.thumb .title a', function(e) {
            e.preventDefault();

            detail.load($(e.target).closest('.thumb'));

            return false;
        });

    });
});


define(['jquery', 'app/api', 'app/ui/controls', 'app/ui/list', 'app/util/sticky'], function ($, api, controls, List, Sticky) {

    var Paper = function (options) {

        this.config = $.extend({
            debug: true,
            element: '#detail',
            document: '.document',
            actions: '.inner > .actions',
            neighbors: '.neighbors-thumbs'
        }, options);

        this.$element = $(this.config.element);
        this.$document = this.$element.find(this.config.document);

        this.id = null;
        this.$actions = null;
        this.actions = null;

        this.neighbors = null;
    };

    function toggleLibraryAddOrTrash($button, added) {
        if (added) {
            $button
                .removeClass('detail-library-add')
                .removeClass('detail-library-restore')
                .addClass('detail-library-trash')
                .find('.eai')
                    .removeClass('eai-library-add')
                    .addClass('eai-library-trash');
        } else {
            $button
                .removeClass('detail-library-trash')
                .addClass('detail-library-restore')
                .find('.eai')
                    .removeClass('eai-library-trash')
                    .addClass('eai-library-add');
        }
    }

    Paper.prototype.log = function() {
        if (this.config.debug) {
            console.log('[Paper] ' + arguments[0], Array.prototype.splice.call(arguments, 1));
        }
    };

    Paper.prototype.getId = function() {
        return this.id;
    };

    Paper.prototype.init = function() {
        this.log('init');

        // API events
        var that = this;
        $('body')
            .on('etalia.publication.pin', function(e, result) {
                if (that.id == result.getId()) {
                    that.$actions
                        .find('.detail-pin')
                        .toggleClass('active', result.isPinned());
                }
            })
            .on('etalia.publication.add', function(e, result) {
                if (that.id == result.getId() && result.isAdded()) {
                    var $button = that.$actions.find('.detail-library-add');
                    toggleLibraryAddOrTrash($button, true);
                }
            })
            .on('etalia.publication.trash', function(e, result) {
                if (that.id == result.getId() && result.isTrashed()) {
                    var $button = that.$actions.find('.detail-library-trash');
                    toggleLibraryAddOrTrash($button, false);
                }
            })
            .on('etalia.publication.restore', function(e, result) {
                if (that.id == result.getId() && result.isAdded()) {
                    var $button = that.$actions.find('.detail-library-restore');
                    toggleLibraryAddOrTrash($button, true);
                }
            });

        this.$document
            .on('click', '.detail-pin', function(e) {
                if (!that.id) throw 'Undefined paper id';

                api.pin(that.id);

                e.stopPropagation();
                return false;
            })
            .on('click', '.detail-ban', function(e) {
                if (!that.id) throw 'Undefined paper id';

                api.ban(that.id);

                e.stopPropagation();
                return false;
            })
            .on('click', '.detail-library-add:visible', function(e) {
                if (!that.id) throw 'Undefined paper id';

                api.add(that.id);

                e.stopPropagation();
                return false;
            })
            .on('click', '.detail-library-trash:visible', function(e) {
                if (!that.id) throw 'Undefined paper id';

                api.trash(that.id);

                e.stopPropagation();
                return false;
            })
            .on('click', '.detail-library-restore:visible', function(e) {
                if (!that.id) throw 'Undefined paper id';

                api.restore(that.id);

                e.stopPropagation();
                return false;
            });

        return this;
    };

    Paper.prototype.update = function() {
        this.log('update');

        this.$actions = this.$element.find(this.config.actions);
        this.id = parseInt(this.$actions.data('paper-id'));

        // Sticky action
        this.actions = new Sticky({
            debug: this.config.debug,
            element: this.$actions,
            parent: this.$document,
            top: 20,
            bottom: 20
        });
        this.actions.enable();

        this.neighbors = new List({
            debug: this.config.debug,
            element: this.config.element + ' ' + this.config.neighbors
        }).init().load(controls.getStates());

        return this;
    };

    Paper.prototype.clear = function() {
        if (this.actions) {
            this.actions.disable();
            this.actions = null;
        }

        if (this.$actions) {
            this.$actions = null;
        }

        if (this.id) {
            this.id = null;
        }
    };

    return Paper;
});
